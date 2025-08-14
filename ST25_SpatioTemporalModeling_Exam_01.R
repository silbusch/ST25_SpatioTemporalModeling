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
library(maxnet)

#---- Download of "Oreamnos americanus (Blainville, 1816)" ---------------------

#**############################################################################*
#** Either first download of "Oreamnos americanus (Blainville, 1816)"...*

# path <- setwd("....")
#
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

#**...or if the data already exists locally: *
path <- "C:/users/Duck/Documents/Studium/EAGLE/2_semester/3_Spatio_Temporal_Modelling/ST25_SpatioTemporalModeling_Exam/ST25_SpatioTemporalModeling/data/0003051-250811113504898.zip"
tmpdir   <- tempdir()

# unpack occurrence.txt
unzip(path, files = "occurrence.txt", exdir = tmpdir, overwrite = TRUE)
# load data
occ_data <- read_tsv(file.path(tmpdir, "occurrence.txt"), show_col_types = FALSE)
#**############################################################################*
#*



#---- Clean and filtering Data -------------------------------------------------

# Keeping only entries with accuracy of 1 km, and deleting entries of preserved specimen
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




#---- Download bioclimatic variables and altitude data -------------------------

path_data <- "C:/users/Duck/Documents/Studium/EAGLE/2_semester/3_Spatio_Temporal_Modelling/ST25_SpatioTemporalModeling_Exam/ST25_SpatioTemporalModeling/data"

#**############################################################################*
#** Either download bioclimatic variables and altitude data data for the first time*
# bioclim_aoi_file <- file.path(path_data, "bioclim_AOI_2p5m.tif")
# # DOwnload bioclimatic variables and altitude data
# bioclim_data <- worldclim_global(var = "bio",
#                                  res = 2.5,
#                                  path = "data/")
# alt <- elevation_global(res = 2.5, path = "data/")
# bioclim_data <- c(bioclim_data, alt)
# 
# # Crop
# bioclim_data <- crop(x = bioclim_data, y = aoi_ext)
# 
# # Save as GeoTIFF
# writeRaster(bioclim_data,
#             filename = bioclim_aoi_file,
#             overwrite = TRUE)

#** Or load data if already downloaded*
bioclim_aoi_file <- file.path(path_data, "bioclim_AOI_2p5m.tif")
bioclim_data <- rast(bioclim_aoi_file)
#**############################################################################*

plot(bioclim_data[[1]])

## Create 5000 pseudo-absence background points
set.seed(20)
background <- spatSample(x = bioclim_data, size = 5000, values = FALSE, na.rm = TRUE, xy = TRUE)

# Map the background points
points(background,
       col = "grey30",
       pch = 1,
       cex = 0.75)

# Merge background and occurence to one data set
goat_pres <- data.frame(occ_data_filtered, occ = 1) # 1= presence
goat_seed <- background %>% as.data.frame() %>%
  rename(
    decimalLongitude = x,
    decimalLatitude = y) %>% 
    mutate(occ = 0)

goat <- rbind(goat_pres,goat_seed)

# extract variables for each point  
e <- extract(x = bioclim_data, y = goat[,c("decimalLongitude","decimalLatitude")], cells = TRUE)
goat <- cbind(goat,e)

# remove points that fall within the same cell, to only have one point per raster-cell (because the values are similar)
goat <-goat[!duplicated(goat$cell),]

# plot observation and random seed points
ggplot() +
  geom_sf(data = world_sf, color = "grey50", fill = "white", linewidth = 0.2) +
  geom_point(data = goat, aes(x = decimalLongitude, y = decimalLatitude, color = factor(occ)),size = 0.5) +
  scale_color_manual(values = c("0" = "grey30", "1" = "red"), labels = c("Random", "Presence")) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  labs(title = "Presences & Random")




#---- Check bioclimatic variables ----------------------------------------------

# Check correlation of our variables
valnum <- lapply(bioclim_data, as.data.frame ) 
valnum <- do.call(cbind,valnum)

# using pearson correlation, since the variables are metric and continuous
cor_matrix <- cor(valnum, use = "complete.obs", method = "pearson")
# plot correlation matrix
corrplot(cor_matrix, method = "number", type = "upper", tl.cex = 0.4, number.cex = 0.7)

# remove correlations >= 0.7, to reduce multicollinearity risk and make the model more stable
to_remove <- findCorrelation(cor_matrix, cutoff = 0.70, names = TRUE)
all_vars <- names(bioclim_data)
print(all_vars)

# remove strong correlated cariables but keep elevation
to_remove_but_keep_elev <- setdiff(to_remove, "wc2.1_2.5m_elev")
vars_with_elev <- setdiff(all_vars, to_remove_but_keep_elev)
print(vars_with_elev)

# remove elevation from to_remove for a model run without elevation data
vars_without_elev <- setdiff(all_vars, union(to_remove, "wc2.1_2.5m_elev"))
print(vars_without_elev)



#---- Generating test- and training data ---------------------------------------

# generate 70% train and 30% test-data 
# generating random training sample
train_data <- sample(seq_len(nrow(goat)), size=round(0.7*nrow(goat)))

goat_train <-goat[train_data,]
goat_test <- goat[-train_data,]

# checking the balance between presence and background points
print(table(goat_train$occ))
print(table(goat_test$occ))




#---- Decision with or without elevation data ----------------------------------

#**############################################################################*
#**Decide whether the following models run with or without elevation**
#** FALSE: without elevation*
#** TRUE: with elevation*

elev <- FALSE  

if (elev) {
  filtered_var <- vars_with_elev
} else {
  filtered_var <- vars_without_elev
}
#**############################################################################*


################################################################################
#-------------------------------------------------------------------------------
#---- Generalized Linear Model -------------------------------------------------
#-------------------------------------------------------------------------------

# generate model on training data
var_glm <- as.formula(paste("occ ~", paste(filtered_var, collapse = " + ")))

# selection of predictors
model_glm <- step(glm(var_glm, family = binomial(link = "logit"), data = goat_train))

summary(model_glm)

# predict based on test data
test_preds <- predict(model_glm, newdata = goat_test, type = "response")
train_preds <- predict(model_glm, newdata = goat_train, type = "response")



#---- Check for overfitting ----------------------------------------------------

# compute AUC and create ROC
roc_test <- roc(goat_test$occ, test_preds)
roc_train <- roc(goat_train$occ, train_preds)

# TODO: Using the area between the curves instead of substraction???
# the gap of the ROC-curve between training and test is an indicator for overfitting: big gab --> overfitting
auc_train_glm <- auc(roc_train)
auc_test_glm  <- auc(roc_test)
print(c(AUC_train=as.numeric(auc_train_glm), AUC_test=as.numeric(auc_test_glm),
        gap=as.numeric(auc_train_glm-auc_test_glm)))

# plot train and test-ROC curve and write down AUC values and gap between ROC-curves
plot(roc_train, col = "blue")
lines(roc_test, col = "red")
legend(
  "bottomright", legend = c(
    paste0("Train (AUC = ", round(auc_train_glm, 3), ")"),
    paste0("Test  (AUC = ", round(auc_test_glm, 3), 
           ", Gap = ", round(gap_glm, 3), ")")
  ),
  col = c("blue", "red"), lwd = 2)

#---- Predict to raster --------------------------------------------------------

# extract relevant predictors
bioclim_sel <- bioclim_data[[filtered_var]]

# using only predictors, that were actually used in the model
pred_in_model <- attr(terms(model_glm), "term.labels")
bioclim_sel   <- bioclim_sel[[intersect(pred_in_model, names(bioclim_sel))]]

# prediction based on GLM-model
model_glm_pred <- terra::predict(object = bioclim_sel,model  = model_glm, type   = "response",na.rm  = TRUE)

plot(model_glm_pred, main = "Predicted Habitat Suitability (GLM)")




################################################################################
#-------------------------------------------------------------------------------
#---- Maximum Entropy Model ----------------------------------------------------
#-------------------------------------------------------------------------------

# filtering for only columns, which were used in the GLM-model and deleting a few columns with NA-values
pred_cols <- intersect(filtered_var, names(goat))
mx_train <- goat_train[, c("occ", pred_cols)]
mx_train <- mx_train[complete.cases(mx_train[, pred_cols]), ]
mx_test  <- goat_test[,  c("occ", pred_cols)]
mx_test  <- mx_test[complete.cases(mx_test[, pred_cols]), ]

#
fmx <- maxnet.formula(p = mx_train$occ, data = mx_train[, pred_cols], classes = "lqph")

# maxnet model 
model_maxnet <- maxnet(p= mx_train$occ, data = mx_train[, pred_cols], f= fmx)

summary(model_maxnet)

pred_train_mx <- predict(model_maxnet, newdata = mx_train[, pred_cols], type = "cloglog")
pred_test_mx  <- predict(model_maxnet, newdata = mx_test[,  pred_cols], type = "cloglog")


#---- Check for overfitting ----------------------------------------------------

# compute AUC and create ROC
roc_train_mx <- pROC::roc(mx_train$occ, as.numeric(pred_train_mx))
roc_test_mx  <- pROC::roc(mx_test$occ,  as.numeric(pred_test_mx))

# TODO: Using the area between the curves instead of substraction???
# the gap of the ROC-curve between training and test is an indicator for overfitting: big gab --> overfitting
auc_train_mx <- pROC::auc(roc_train_mx)
auc_test_mx  <- pROC::auc(roc_test_mx)
gap_mx       <- as.numeric(auc_train_mx - auc_test_mx)

# plot train and test-ROC curve and write down AUC values and gap between ROC-curves
plot(roc_train_mx, col = "blue")
lines(roc_test_mx, col = "red")
legend(
  "bottomright", 
  legend = c(
    paste0("Train (AUC = ", round(auc_train_mx, 3), ")"),
    paste0("Test  (AUC = ", round(auc_test_mx, 3),")"),
    paste0("Gap   = ", round(gap_mx, 3))),
  col = c("blue", "red", NA),
  lwd = c(2, 2, NA),
  bty = "n")


#---- Predict to raster --------------------------------------------------------

# extract relevant predictors
bioclim_mx <- bioclim_data[[pred_cols]]

# MaxEnt model using for prediction for each raster zell (Values between 0 and 1 (cloglog))
model_maxnet_pred <- terra::predict(bioclim_mx, model_maxnet, type = "cloglog", na.rm = TRUE)
plot(model_maxnet_pred, main = "Predicted Habitat Suitability (MaxNet)")




################################################################################
#-------------------------------------------------------------------------------
#---- Random Forest Model ------------------------------------------------------
#-------------------------------------------------------------------------------










par(mfrow = c(2, 4))
plot(model_glm_pred, main = "GLM Habitat Suitability")
plot(model_maxnet_pred, main = "MaxEnt Habitat Suitability")
plot(..., main = "...")
plot(..., main = "...")
par(mfrow = c(1, 1))
