################################################################################
#### Habitat modeling - Oreamnos americanus ####################################
################################################################################

library(readr)
library(rgbif)
library(dplyr)
library(maps)
library(terra)
library(geodata)
library(sf)
library(ggplot2)
library(lattice)
library(corrplot)
library(caret)
library(pROC)
library(predicts)

#---- Download data of "Oreamnos americanus (Blainville, 1816)" ----------------

path <- setwd("....")

# # Download Data
# key <- name_backbone(name = "Oreamnos americanus")$usageKey
# 
# download_key <- occ_download(
#   pred("taxonKey", key),
#   pred("hasCoordinate", TRUE),
#   user = "...",
#   pwd = "...",
#   email = "..."
# )
# 
# dl <- occ_download_get(download_key, overwrite = TRUE)
# 
# occ_data <- occ_download_import(dl)

# If data is already local
path <- "C:/users/Duck/Documents/Studium/EAGLE/2_semester/3_Spatio_Temporal_Modelling/ST25_SpatioTemporalModeling_Exam/ST25_SpatioTemporalModeling/data/0003051-250811113504898.zip"
tmpdir   <- tempdir()

# occurrence.txt entpacken
unzip(path, files = "occurrence.txt", exdir = tmpdir, overwrite = TRUE)

# Einlesen
occ_data <- read_tsv(file.path(tmpdir, "occurrence.txt"), show_col_types = FALSE)

#---- Clean and filtering Data -------------------------------------------------
# Keeping only entries with accuracy of 1 km,
# and deleting entries of preserved specimen
occ_data_filtered <- occ_data %>%
  filter(
    (is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= 1000) &
      !basisOfRecord %in% c("PRESERVED_SPECIMEN")&
      decimalLongitude < -62)

# remove unnecessary columns
occ_data_filtered <- occ_data_filtered %>% select(c("decimalLongitude","decimalLatitude"))

#---- Define AOI ---------------------------------------------------------------
# Extent of AOI
xmin <- -160
xmax <-  -100
ymin <-   35
ymax <-   65

aoi_ext <- ext(xmin, xmax, ymin, ymax)

# Preparing world map for AOI plot
world <- world(resolution = 3, path = "data/")
world_crop <- crop(world, aoi_ext)
world_sf <- st_as_sf(world_crop)

# Plot AOI
ggplot() +
  geom_sf(data = world_sf, color = "grey50", linewidth = 0.2) +
  geom_point(data = occ_data_filtered,
             aes(x = decimalLongitude, y = decimalLatitude),
             size = 0.5, color = "red") +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  labs (title = "Area of interest")

#---- Bioclimatic variables, altitude data -------------------------------------
path_data <- "C:/users/Duck/Documents/Studium/EAGLE/2_semester/3_Spatio_Temporal_Modelling/ST25_SpatioTemporalModeling_Exam/ST25_SpatioTemporalModeling/data"

# Bioclimatic variables
#TODO: Download 0.5 --> 1km
# Download for USA and Canada

# bio_us  <- worldclim_country("USA", var="bio",  res=0.5, path=path_data)  |> crop(geographic_extent)
# bio_can <- worldclim_country("CAN", var="bio",  res=0.5, path=path_data)  |> crop(geographic_extent)
# bio     <- mosaic(bio_us, bio_can, fun="mean")

# bio_us <- worldclim_country(country = "USA", var = "bio", res = 2.5, path = path_data)
# bio_ca <- worldclim_country(country = "CAN", var = "bio", res = 2.5, path = path_data)
# 
# bio    <- mosaic(bio_us, bio_ca, fun = "mean")
# bio    <- crop(bio, aoi_ext)
# 
# # Altitude data
# elev_us <- worldclim_country(country = "USA", var = "elev", res = 2.5, path = path_data)
# elev_ca <- worldclim_country(country = "CAN", var = "elev", res = 2.5, path = path_data)
# elev    <- mosaic(elev_us, elev_ca, fun = "mean")
# elev    <- crop(elev, aoi_ext)

################################################################################
# DOwnload bioclimatic variables and altitude data
bioclim_data <- worldclim_global(var = "bio",
                                 res = 2.5,
                                 path = "data/")
alt <- elevation_global(res = 2.5, path = "data/")
bioclim_data <- c(bioclim_data, alt)

bioclim_data <- crop(x = bioclim_data, y = aoi_ext)
plot(bioclim_data[[1]])


## Create background points
set.seed(20)

background <- spatSample(x = bioclim_data,
                         size = 5000,    # generate 5000 pseudo-absence points
                         values = FALSE, # don't need values
                         na.rm = TRUE,   
                         xy = TRUE)

# Map the background points
points(background,
       col = "grey30",
       pch = 1,
       cex = 0.75)

## Merge background and occurence to one data set
goat <- data.frame(occ_data_filtered, occ = 1) # 1= presence
goat <- background %>% as.data.frame() %>%
  rename(
    decimalLongitude = x,
    decimalLatitude = y
  ) %>% mutate(
    occ = 0
  ) %>% 
  rbind(goat)


e <- extract(x = bioclim_data, y = goat[,c("decimalLongitude","decimalLatitude")], cells = TRUE)
goat <- cbind(goat, e)
## remove points that fall within the same cell, to only have one point per raster-cell (becuase the values are similar)
goat <-goat[!duplicated(goat$cell),]

# Check correlation of our variables
valnum <- lapply(bioclim_data, as.data.frame ) 
valnum <- do.call(cbind,valnum)

cor_matrix <- cor(valnum, use = "complete.obs", method = "pearson")
corrplot(cor_matrix, method = "number", type = "upper", tl.cex = 0.4, number.cex = 0.7)

to_remove <- findCorrelation(cor_matrix, cutoff = 0.75, names = TRUE)

# remove elevation from to_remove because it makes a difference (try once without and once with elevation in)
to_remove <- setdiff(to_remove, "wc2.1_2.5m_elev")
bioclim_data_pred <- bioclim_data[[!names(bioclim_data) %in% to_remove]]

#########
###### Make a model with training and evaluation data
#########

## Model evaluation and validation
# Generate indexes for generating training and test/validation data - 70% of the data should be train and 30 for testing. You can also use other thresholds
train_data <- sample(seq_len(nrow(goat)), size=round(0.7*nrow(goat)))

goat_train <-goat[train_data,]
goat_test <- goat[-train_data,]

# generate model on training data
model_glm = step(glm(occ ~ wc2.1_2.5m_bio_2 + wc2.1_2.5m_bio_3 + wc2.1_2.5m_bio_6 + wc2.1_2.5m_bio_8 + 
                       wc2.1_2.5m_bio_9 + wc2.1_2.5m_bio_13 +wc2.1_2.5m_bio_15 + wc2.1_2.5m_bio_18 + wc2.1_2.5m_bio_19 + wc2.1_2.5m_elev, family=binomial(link=logit), data= goat_train))
summary(model_glm)

# predict based on test data
test_preds <- predict(model_glm, newdata = goat_test, type = "response")

# Plot ROC curve and compute AUC
roc_obj <- roc(goat_test$occ, test_preds)
plot(roc_obj, main = paste("AUC =", round(auc(roc_obj), 3)))
