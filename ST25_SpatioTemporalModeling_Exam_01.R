################################################################################
#### Habitat modeling - Oreamnos americanus ####################################
################################################################################

library(readr)
library(rgbif)
library(dplyr)
library(terra)
library(geodata)
library(sf)
library(ggplot2)
library(lattice)
library(corrplot)
library(caret)
library(pROC)
library(maxnet)
library(randomForest)
library(mgcv)
library(splines)
library(gam)
library(openxlsx)
library(remotes)
library(scales)
remotes::install_github("rvalavi/blockCV", dependencies = TRUE)
library(blockCV)

#---- Download of "Oreamnos americanus (Blainville, 1816)" ---------------------

#*******************************************************************************
#** Either first download of "Oreamnos americanus (Blainville, 1816)"...*

path <- setwd("...")

# Download Data
key <- name_backbone(name = "Oreamnos americanus")$usageKey

download_key <- occ_download(
  pred("taxonKey", key),
  pred("hasCoordinate", TRUE),
  user = "...",
  pwd = "...",
  email = "..."
)

dl <- occ_download_get(download_key, overwrite = TRUE)

occ_data <- occ_download_import(dl)

# *******************************************************************************
# #**...or if the data already exists locally: *
# 
# path <- "C:/users/Duck/Documents/Studium/EAGLE/2_semester/3_Spatio_Temporal_Modelling/ST25_SpatioTemporalModeling_Exam/ST25_SpatioTemporalModeling/data/0003051-250811113504898.zip"
# 
# path <- "C:/.../0003051-250811113504898.zip"
# tmpdir   <- tempdir()
# 
# # unpack occurrence.txt
# unzip(path, files = "occurrence.txt", exdir = tmpdir, overwrite = TRUE)
# # load data
# occ_data <- read_tsv(file.path(tmpdir, "occurrence.txt"), show_col_types = FALSE)
# 
# *******************************************************************************




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
xmin <- -154
xmax <-  -103
ymin <-   37
ymax <-   63

aoi_ext <- ext(xmin, xmax, ymin, ymax)

# Preparing world map for AOI plot
world <- world(resolution = 3, path = "data/")
world_crop <- crop(world, aoi_ext)
world_sf <- st_as_sf(world_crop)

# Plot AOI
ggplot() +
  geom_sf(data = world_sf, color = "grey50", linewidth = 0.2) +
  geom_point(data = occ_data_filtered,
             aes(x = decimalLongitude, y = decimalLatitude, color = "Observed occurrence"),
             size = 0.5) +
  scale_color_manual(name = "Legend", values = c("Observed occurrence" = "red")) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  labs (title = "Area of interest of the habitat of Oreamnos americanus",
        x= "Longitude",
        y= "Latitude")




#---- Download bioclimatic variables and altitude data -------------------------

#*******************************************************************************
#** Either download bioclimatic variables and altitude data for the first time*

# path_data <- "C:/users/Duck/Documents/Studium/EAGLE/2_semester/3_Spatio_Temporal_Modelling/ST25_SpatioTemporalModeling_Exam/ST25_SpatioTemporalModeling/data"

# location where the bioclimatic- and elevation data should be stored
path_data <- "C:/...."

bioclim_aoi_file <- file.path(path_data, "bioclim_AOI_2p5m.tif")

# DOwnload bioclimatic variables and altitude data
bioclim_data <- worldclim_global(var = "bio",
                                 res = 2.5,
                                 path = "data/")

alt <- elevation_global(res = 2.5, path = "data/")

bioclim_data <- c(bioclim_data, alt)

# Crop
bioclim_data <- crop(x = bioclim_data, y = aoi_ext)

# Save as GeoTIFF
writeRaster(bioclim_data,
            filename = bioclim_aoi_file,
            overwrite = TRUE)

#** Or load data if already downloaded*
#*
# bioclim_aoi_file <- file.path(path_data, "bioclim_AOI_2p5m.tif")
# bioclim_data <- rast(bioclim_aoi_file)

#*******************************************************************************

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
table(goat$occ)

# plot observation and random seed points
ggplot() +
  geom_sf(data = world_sf, color = "grey50", fill = "white", linewidth = 0.2) +
  geom_point(data = goat, aes(x = decimalLongitude, y = decimalLatitude, color = factor(occ)),size = 0.5) +
  scale_color_manual(name= "Legend", values = c("0" = "grey30", "1" = "red"), labels = c("Random", "Observed Occurrence")) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  labs(title = "Observed occurrence & Random",
       x= "Longitude",
       y= "Latitude")




#---- Check bioclimatic variables ----------------------------------------------

# Check correlation of our variables
valnum <- lapply(bioclim_data, as.data.frame )
valnum <- do.call(cbind,valnum)

# using pearson correlation, since the variables are metric and continuous
cor_matrix <- cor(valnum, use = "complete.obs", method = "pearson")
# plot correlation matrix
corrplot(cor_matrix, method = "number", type = "upper", tl.cex = 0.4, number.cex = 0.7)

# remove correlations >= 0.70, to reduce multicollinearity risk and make the model more stable
to_remove <- findCorrelation(cor_matrix, cutoff = 0.7, names = TRUE)
all_vars <- names(bioclim_data)
print(all_vars)

# remove strong correlated cariables but keep elevation
to_remove_but_keep_elev <- setdiff(to_remove, "wc2.1_2.5m_elev")
vars_with_elev <- setdiff(all_vars, to_remove_but_keep_elev)
print(vars_with_elev)

# remove elevation from to_remove for a model run without elevation data
vars_without_elev <- setdiff(all_vars, union(to_remove, "wc2.1_2.5m_elev"))
print(vars_without_elev)


# example variable selection option above leaded to AUC for Train and Test between 0.9 and 0.99 for all models 
# all_vars <- names(bioclim_data)
# 
# keep_vars <- c(
#   "wc2.1_2.5m_bio_1",
#   "wc2.1_2.5m_bio_2",
#   "wc2.1_2.5m_bio_3",
#   "wc2.1_2.5m_bio_8",
#   "wc2.1_2.5m_bio_9",
#   "wc2.1_2.5m_bio_10",
#   "wc2.1_2.5m_bio_11",
#   "wc2.1_2.5m_bio_13",
#   "wc2.1_2.5m_bio_14",
#   "wc2.1_2.5m_bio_18",
#   "wc2.1_2.5m_bio_19")
# 
# keep_vars <- all_vars[all_vars %in% keep_vars]
# elev_name <- "wc2.1_2.5m_elev"
# 
# vars_with_elev    <- unique(c(keep_vars, intersect(elev_name, all_vars)))
# 
# # remove elevation from to_remove for a model run without elevation data
# vars_without_elev <- setdiff(keep_vars, elev_name)


#---- Decision with or without elevation data ----------------------------------

#*******************************************************************************
#** Decide whether the following models run with or without elevation *
#** FALSE: without elevation*
#** TRUE: with elevation*

elev <- FALSE 

if (elev) {
  filtered_var <- vars_with_elev
} else {
  filtered_var <- vars_without_elev
}
#*******************************************************************************

#---- Generating test- and training data ---------------------------------------

#** The 70/30 split has the risk of spatial autocorrelation, so blockCV is used*

# # generate 70% train and 30% test-data 
# # generating random training sample
# train_data <- sample(seq_len(nrow(goat)), size=round(0.7*nrow(goat)))
# 
# goat_train <-goat[train_data,]
# goat_test <- goat[-train_data,]
# 
# # checking the balance between presence and background points
# print(table(goat_train$occ))
# print(table(goat_test$occ))


#---- Block Cross Validation ---------------------------------------------------

# To account for the effect of spatial autocorrelation due to closely located
# test and training data, blockCV is used to create spatially seperated folds. 

goat$decimalLongitude <- as.numeric(goat$decimalLongitude)
goat$decimalLatitude  <- as.numeric(goat$decimalLatitude)

# points as sf-object and taking crs
goat_sf <- st_as_sf(
  goat,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326,
  remove = FALSE)

bioclim_sel <- bioclim_data[[filtered_var]]

# setting block size manually, cause automatic setting use of the "cv_spatial_autocor"
# function led to overly large blocks, also after repeatedly trying out different variables
size_m <- 200000  

# The presence and 5000 pseudo points are not distributed randomly, but are now 
# divided up based on the hexagons in training or test data.
old_seed <- .Random.seed
set.seed(50)
cv <- blockCV::cv_spatial(
  x = goat_sf,
  column = "occ",
  r = bioclim_sel,
  k = 10,
  size = size_m,
  hexagon = TRUE,
  selection = "random",
  iteration = 100,
  plot = TRUE
)
.Random.seed <- old_seed

folds <- cv$folds_list

################################################################################
##         Model Area                                                         ##
################################################################################

# result dataframe for all folds
results_all <- data.frame()

# create list to save prediction for every fold prediction of every model
mean_preds <- list(GLM=NULL, GAM=NULL, MaxNet=NULL, RF=NULL, Ensemble=NULL)

r_all <- bioclim_data[[filtered_var]]

#** !!! the loop may take several minutes to complete !!!*

for (i in seq_along(folds)) {
  # outputs the current fold to show progress
  cat("\n--- Fold", i, "---\n")
  
#---- create train and test data -----------------------------------------------
  train_idx <- unlist(folds[[i]][[1]])
  test_idx  <- unlist(folds[[i]][[2]])
  
  goat_train <- goat[train_idx, ]
  goat_test  <- goat[test_idx, ]

################################################################################
#-------------------------------------------------------------------------------
#---- Generalized Linear Model -------------------------------------------------
#------------------------------------------------------------------------------- 
  
  # formular for GLM
  var_glm <- as.formula(paste("occ ~", paste(filtered_var, collapse = " + ")))
  # stepwise selektion of variables
  model_glm <- step(glm(var_glm, family = binomial(link = "logit"), data = goat_train))
  
  # prediction for train and test data
  pred_train_glm <- predict(model_glm, newdata = goat_train, type = "response")
  pred_test_glm  <- predict(model_glm, newdata = goat_test, type = "response")
  
#---- Check for overfitting ----------------------------------------------------
  # compute AUC and create ROC
  auc_train_glm <- pROC::auc(goat_train$occ, pred_train_glm)
  auc_test_glm  <- pROC::auc(goat_test$occ,  pred_test_glm)
  # the gap of the ROC-curve between training and test is an indicator for overfitting: big gab --> overfitting
  gap_glm <- as.numeric(auc_train_glm - auc_test_glm)
  
#---- Predict to raster --------------------------------------------------------  
  # using only predictors, that were actually used in the model
  pred_in_model <- attr(terms(model_glm), "term.labels")
  r_glm <- r_all[[intersect(pred_in_model, names(r_all))]]
  # prediction to raster based on GLM-model
  pred_r_glm <- terra::predict(r_glm, model_glm, type = "response", na.rm = TRUE)
  
  
  

################################################################################
#-------------------------------------------------------------------------------
#---- Generalized Additive Model -----------------------------------------------
#-------------------------------------------------------------------------------
  
  # formular for GAM
  var_gam <- as.formula(paste("occ ~", paste0("s(", filtered_var, ")", collapse = " + ")))
  model_gam <- mgcv::gam(var_gam, data = goat_train, family = binomial)
  
  # prediction for train and test data
  pred_train_gam <- predict(model_gam, newdata = goat_train, type = "response")
  pred_test_gam  <- predict(model_gam, newdata = goat_test, type = "response")
#---- Check for overfitting ----------------------------------------------------
  
  auc_train_gam <- pROC::auc(goat_train$occ, pred_train_gam)
  auc_test_gam  <- pROC::auc(goat_test$occ,  pred_test_gam)
  gap_gam <- as.numeric(auc_train_gam - auc_test_gam)

#---- Predict to raster --------------------------------------------------------
  vars_in_gam <- setdiff(all.vars(var_gam), "occ")
  r_gam <- r_all[[intersect(vars_in_gam, names(r_all))]]
  pred_r_gam <- terra::predict(r_gam, model_gam, type = "response", na.rm = TRUE)


  
  
################################################################################
#-------------------------------------------------------------------------------
#---- Maximum Entropy Model ----------------------------------------------------
#-------------------------------------------------------------------------------
  
  # filtering for only columns, which were used in the GLM-model and deleting a few columns with NA-values
  pred_cols <- intersect(filtered_var, names(goat_train))
  
  # make again sure that NAs are removed, otherwise the model wil fail
  mx_train <- goat_train[, c("occ", pred_cols)]
  mx_train <- mx_train[complete.cases(mx_train), ]
  mx_test  <- goat_test[,  c("occ", pred_cols)]
  mx_test  <- mx_test[complete.cases(mx_test), ]
  
  # formular for maxnet model
  fmx <- maxnet.formula(p = mx_train$occ, data = mx_train[, pred_cols], classes = "lqph")
  # maxnet model 
  model_maxnet <- maxnet(p= mx_train$occ, data = mx_train[, pred_cols], f= fmx)
  
  pred_train_mx <- predict(model_maxnet, newdata = mx_train[, pred_cols], type = "cloglog")
  pred_test_mx  <- predict(model_maxnet, newdata = mx_test[,  pred_cols], type = "cloglog")
  
#---- Check for overfitting ----------------------------------------------------
  
  auc_train_mx <- pROC::auc(mx_train$occ, as.numeric(pred_train_mx))
  auc_test_mx  <- pROC::auc(mx_test$occ,  as.numeric(pred_test_mx))
  gap_mx <- as.numeric(auc_train_mx - auc_test_mx)
  
  # extract relevant predictors
  r_mx <- r_all[[pred_cols]]
  
  # MaxEnt model using for prediction for each raster zell (Values between 0 and 1 (cloglog))
  pred_r_mx <- terra::predict(r_mx, model_maxnet, type = "cloglog", na.rm = TRUE)

  
  
 
################################################################################
#-------------------------------------------------------------------------------
#---- Random Forest Model ------------------------------------------------------
#-------------------------------------------------------------------------------
# Random Forest model with Bootstrap-Sampling 
  rf_train <- goat_train[, c("occ", pred_cols)]
  rf_test  <- goat_test[,  c("occ", pred_cols)]
  rf_train <- rf_train[complete.cases(rf_train), ]
  rf_test  <- rf_test[complete.cases(rf_test), ]
  rf_train$occ <- factor(rf_train$occ, levels = c(0, 1))
  rf_test$occ  <- factor(rf_test$occ,  levels = c(0, 1))
  
  #RF- model with 500 trees (more could lead to overfitting), bootstrap per tree
  model_rf <- randomForest(occ ~ ., data = rf_train, ntree = 500, importance = FALSE)
  
#---- Check for overfitting ----------------------------------------------------
  
  pred_train_rf <- predict(model_rf, type = "prob")[, "1"]
  pred_test_rf  <- predict(model_rf, newdata = rf_test, type = "prob")[, "1"]
  auc_train_rf <- pROC::auc(rf_train$occ, as.numeric(pred_train_rf))
  auc_test_rf  <- pROC::auc(rf_test$occ,  as.numeric(pred_test_rf))
  gap_rf <- as.numeric(auc_train_rf - auc_test_rf)

#---- Predict to raster --------------------------------------------------------
  
  r_rf <- r_all[[pred_cols]]
  
  pred_r_rf <- terra::predict(r_rf, model_rf, type = "prob", index = 2, na.rm = TRUE)
  
  
  
  
################################################################################
#-------------------------------------------------------------------------------
#---- Ensemble Model -----------------------------------------------------------
#-------------------------------------------------------------------------------
  
  # mean value for all predicted test and train data
  ens_train <- rowMeans(cbind(
    as.numeric(pred_train_glm),
    as.numeric(pred_train_gam),
    as.numeric(pred_train_mx),
    as.numeric(pred_train_rf)), na.rm = TRUE)
  
  ens_test <- rowMeans(cbind(
    as.numeric(pred_test_glm),
    as.numeric(pred_test_gam),
    as.numeric(pred_test_mx),
    as.numeric(pred_test_rf)), na.rm = TRUE)
  
  # Roc, Auc and Gap
  auc_train_ens <- pROC::auc(goat_train$occ, ens_train)
  auc_test_ens  <- pROC::auc(goat_test$occ,  ens_test)
  gap_ens <- as.numeric(auc_train_ens - auc_test_ens)
  
  #All grids of the model predictions for the [i]-fold run are added together 
  # and their cell values are divided by four to average the result.
  pred_r_ens <- (pred_r_glm + pred_r_gam + pred_r_mx + pred_r_rf) / 4
  
  
  
  
################################################################################
#-------------------------------------------------------------------------------
#---- Saving all results -------------------------------------------------------
#------------------------------------------------------------------------------- 

  # saving all AUC and Gap results 
  results_all <- rbind(results_all,
    data.frame(Fold = i, Model = "GLM", AUC_train = auc_train_glm, AUC_test = auc_test_glm, Gap = gap_glm),
    data.frame(Fold = i, Model = "GAM", AUC_train = auc_train_gam, AUC_test = auc_test_gam, Gap = gap_gam),
    data.frame(Fold = i, Model = "MaxNet", AUC_train = auc_train_mx,  AUC_test = auc_test_mx,  Gap = gap_mx),
    data.frame(Fold = i, Model = "RF", AUC_train = auc_train_rf,  AUC_test = auc_test_rf,  Gap = gap_rf),
    data.frame(Fold = i, Model = "Ensemble", AUC_train = auc_train_ens, AUC_test = auc_test_ens, Gap = gap_ens))
 
  # saving alls [i] rasters, add them together
  if (is.null(mean_preds$GLM)) {
    mean_preds$GLM      <- pred_r_glm
    mean_preds$GAM      <- pred_r_gam
    mean_preds$MaxNet   <- pred_r_mx
    mean_preds$RF       <- pred_r_rf
    mean_preds$Ensemble <- pred_r_ens
  } else {
    mean_preds$GLM      <- mean_preds$GLM      + pred_r_glm
    mean_preds$GAM      <- mean_preds$GAM      + pred_r_gam
    mean_preds$MaxNet   <- mean_preds$MaxNet   + pred_r_mx
    mean_preds$RF       <- mean_preds$RF       + pred_r_rf
    mean_preds$Ensemble <- mean_preds$Ensemble + pred_r_ens
  }
}


################################################################################
##         Table with all AUC Values                                          ##
################################################################################

summary_results <- results_all %>%
  group_by(Model) %>%
  summarise(
    Mean_AUC_train = mean(AUC_train),
    SD_AUC_train   = sd(AUC_train),
    Mean_AUC_test  = mean(AUC_test),
    SD_AUC_test    = sd(AUC_test),
    Mean_Gap       = mean(Gap),
    SD_Gap         = sd(Gap),
    .groups = "drop"  )

print(summary_results)

summary_results_rounded <- summary_results %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

# creating table depending if elev = TRUE or FALSE
file_suffix <- if (elev) "with_elev" else "without_elev"
# individual naming
file_name   <- paste0("AUC_results_", file_suffix, ".xlsx")

# create .xlsx table in working directory
openxlsx::write.xlsx(summary_results_rounded, file = file_name, rowNames = FALSE)
cat("Table got saved as:", file_name, "\n")

################################################################################
##         All predicted rasters                                              ##
################################################################################

# calculating the average prediction raster for each model
k_folds <- length(folds)

for (m in names(mean_preds)) {
  mean_preds[[m]] <- mean_preds[[m]] / k_folds}

plot_suffix <- if (elev) "with_elev" else "without_elev"
mean_plot_file <- paste0("Suitability_models_mean_", plot_suffix, ".png")

# convert raster in dataframe
df_glm <- as.data.frame(mean_preds$GLM, xy = TRUE) %>% mutate(Model = "GLM")
df_gam <- as.data.frame(mean_preds$GAM, xy = TRUE) %>% mutate(Model = "GAM")
df_mx  <- as.data.frame(mean_preds$MaxNet, xy = TRUE) %>% mutate(Model = "MaxNet")
df_rf  <- as.data.frame(mean_preds$RF, xy = TRUE) %>% mutate(Model = "RF")
df_ens <- as.data.frame(mean_preds$Ensemble, xy = TRUE) %>% mutate(Model = "Ensemble")

# rename to get a value column for all data frames
names(df_glm)[3] <- "Suitability"
names(df_gam)[3] <- "Suitability"
names(df_mx)[3]  <- "Suitability"
names(df_rf)[3]  <- "Suitability"
names(df_ens)[3] <- "Suitability"

# combine
df_all <- bind_rows(df_glm, df_gam, df_mx, df_rf, df_ens)

mean_plot <- ggplot(df_all, aes(x = x, y = y, fill = Suitability)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Suitability\nas a habitat", option = "plasma",limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  facet_wrap(~ Model, ncol = 3) +
  coord_equal() +
  labs(title = "Mean Model Predictions of GLM, GAM, MaxNet, RF, and Ensemble",
       x = "Longitude", y = "Latitude")

ggsave(mean_plot_file, plot = mean_plot, width = 12, height = 5, dpi = 300)
cat("Plot got saved as:", mean_plot_file, "\n")




################################################################################
#-------------------------------------------------------------------------------
#---- Threshold Model ----------------------------------------------------------
#-------------------------------------------------------------------------------
# Representation of the habitat when all grid cells with a certain threshold are
# combined from all models. So, the locations where all models predict a very
# high probability of occurrence of the Mountain Goat.
# The condition for identifying a grid cell as a potential habitat is that at
# least three of the four models have a prediction value of >= 0.7.

# threshold
thr <- 0.70

# binary raster
rf_bin <- terra::ifel(mean_preds$RF >= thr, 1, 0)
glm_bin <- terra::ifel(mean_preds$GLM >= thr, 1, 0)
gam_bin <- terra::ifel(mean_preds$GAM >= thr, 1, 0)
mx_bin  <- terra::ifel(mean_preds$MaxNet >= thr, 1, 0)
ens_bin <- terra::ifel(mean_preds$Ensemble >= thr, 1, 0)

# count cells with valid value
hits <- glm_bin + gam_bin + rf_bin + mx_bin
consensus_bin <- terra::ifel(hits >= 3, 1, 0)

# individual suffix
plot_suffix <- if (elev) "with_elev" else "without_elev"
threshold_plot_file <- paste0("Threshold_models_mean_", plot_suffix, ".png")

png(filename = threshold_plot_file, width = 2000, height = 1500, res = 150)

par(mfrow = c(2, 3), oma = c(0, 0, 4, 12))

par(mfrow = c(2, 3))
plot(glm_bin, main = paste0("GLM ≥ ", thr), col = c("grey80","darkgreen"), legend = FALSE)
plot(gam_bin, main = paste0("GAM ≥ ", thr), col = c("grey80","darkgreen"), legend = FALSE)
plot(mx_bin,  main = paste0("MaxNet ≥ ", thr), col = c("grey80","darkgreen"), legend = FALSE)
plot(rf_bin,  main = paste0("RF ≥ ", thr), col = c("grey80","darkgreen"), legend = FALSE)
plot(ens_bin, main = paste0("Ensemble ≥ ", thr), col = c("grey80","darkgreen"), legend = FALSE)
plot(consensus_bin, main = paste0("Consensus ≥ ", thr, " (≥3 models)"), col = c("grey80","red"), legend = FALSE)
par(xpd = NA)

cols <- c("grey80", "red")

legend("right", legend = c("Probably no habitat", "Very likely habitat"), fill = cols,
       border = "black", inset = -0.4, bty = "n", cex = 1)

mtext("Threshold (prediction values ≥ 0.7", outer = TRUE, line = -2, cex = 1.5)
par(mfrow = c(1, 1))
dev.off()

cat("Threshold maps saved as:", threshold_plot_file, "\n")


################################################################################
##         model differences                                                  ##
################################################################################
# comparing the predicted results 

# stacking mean-models
r_stack <- c(mean_preds$GLM, mean_preds$GAM, mean_preds$MaxNet, mean_preds$RF)
names(r_stack) <- c("GLM","GAM","MaxNet","RF")


#---- Standard deviation -------------------------------------------------------
# standard deviation raster for each cell
r_sd <- app(r_stack, fun = sd, na.rm = TRUE)

plot_suffix <- if (elev) "with_elev" else "without_elev"
sd_plot_file <- paste0("SD_between_models_", plot_suffix, ".png")

df_sd <- as.data.frame(r_sd, xy = TRUE, na.rm = TRUE)
names(df_sd)[3] <- "SD"

p_sd <-ggplot(df_sd, aes(x = x, y = y, fill = SD)) +
  geom_raster() +
  scale_fill_viridis_c(
    name   = "Standard\nDeviation",
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, 0.1)
  ) +
  coord_fixed() +
  labs(
    title = "Standard deviation of GLM, GAM, MaxNet, RF",
    subtitle = "Comparison of spatial predictions",
    x = "Longitude", y = "Latitude")

# save in working directory
ggsave(sd_plot_file, plot = p_sd, width = 12, height = 7, dpi = 300)
cat("Standard deviation raster saved as:", sd_plot_file, "\n")

# histogramm
vals <- as.vector(terra::values(r_sd, na.rm = TRUE))

hist_plot_file <- paste0("SD_histogram_", plot_suffix, ".png")

sd_hist <- ggplot(data.frame(SD = vals), aes(x = SD)) +
  geom_histogram(binwidth = 0.01, fill = "grey", color = "black") +
  scale_x_continuous(
    breaks = seq(0, 0.5, by = 0.05),) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Histogram of the standard deviation",
       x = "Standard deviation",
       y = "Number of grid cells")

# save in working directory
png(filename = hist_plot_file, width = 1600, height = 1100, res = 150)
print(sd_hist)
dev.off()

cat("SD histogram saved as:", hist_plot_file, "\n")

#---- spatial correlation ------------------------------------------------------

set.seed(50)
smp_vals <- spatSample(r_stack, size = 50000, method = "random", na.rm = TRUE, as.points = FALSE, values = TRUE)
cmat <- cor(as.data.frame(smp_vals), use = "pairwise.complete.obs")

# save correlation matrix as png
corr_png_file <- paste0("Model_prediction_correlation_", plot_suffix, ".png")

png(filename = corr_png_file, width = 1400, height = 1200, res = 150)

par(oma = c(0, 0, 6, 0), mar = c(3, 3, 1, 6))
corrplot::corrplot(cmat, method = "color", type = "full", addCoef.col = "black",
                   tl.col = "black", tl.srt = 45, number.cex = 0.8, cl.pos = "r")
mtext("Correlation of mean model predictions", side = 3, outer = TRUE, line = 2, cex = 1.4)

dev.off()

cat("Correlation matrix:", corr_png_file, "\n")

#