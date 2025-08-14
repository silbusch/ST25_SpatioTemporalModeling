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
library(randomForest)
library(mgcv)
library(gbm)
library(splines)
library(foreach)
library(gam)
library(openxlsx)

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
             aes(x = decimalLongitude, y = decimalLatitude),
             size = 0.5, color = "red") +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  labs (title = "Area of interest")




#---- Download bioclimatic variables and altitude data -------------------------

path_data <- "C:/users/Duck/Documents/Studium/EAGLE/2_semester/3_Spatio_Temporal_Modelling/ST25_SpatioTemporalModeling_Exam/ST25_SpatioTemporalModeling/data"

#*******************************************************************************
#*                                                                      ********
#*                                                                      ********
#** Either download bioclimatic variables and altitude data for the first time*
# bioclim_aoi_file <- file.path(path_data, "bioclim_AOI_2p5m.tif")
# # DOwnload bioclimatic variables and altitude data
# bioclim_data <- worldclim_global(var = "bio",
#                                  res = 2.5,
#                                  path = "data/")
# alt <- elevation_global(res = 2.5, path = "data/")
# bioclim_data <- c(bioclim_data, alt)
# 
# Crop
# bioclim_data <- crop(x = bioclim_data, y = aoi_ext)
# 
# # Save as GeoTIFF
# writeRaster(bioclim_data,
#             filename = bioclim_aoi_file,
#             overwrite = TRUE)

#** Or load data if already downloaded*
bioclim_aoi_file <- file.path(path_data, "bioclim_AOI_2p5m.tif")
bioclim_data <- rast(bioclim_aoi_file)
#*                                                                      ********
#*                                                                      ********
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

# remove correlations >= 0.75, to reduce multicollinearity risk and make the model more stable
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

#Variable selection option above leaded to AUC for Train and Test between 0.9 and 0.99 for all models 

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

#*******************************************************************************
#*                                                                      ********
#** Decide whether the following models run with or without elevation   ********
#** FALSE: without elevation*
#** TRUE: with elevation*

elev <- FALSE  

if (elev) {
  filtered_var <- vars_with_elev
} else {
  filtered_var <- vars_without_elev
}
#*                                                                      ********
#*******************************************************************************

################################################################################
##         Model Area                                                         ##
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

# the gap of the ROC-curve between training and test is an indicator for overfitting: big gab --> overfitting
auc_train_glm <- auc(roc_train)
auc_test_glm  <- auc(roc_test)
print(c(AUC_train=as.numeric(auc_train_glm), AUC_test=as.numeric(auc_test_glm),
        gap=as.numeric(auc_train_glm-auc_test_glm)))
gap_glm <- as.numeric(auc_train_glm - auc_test_glm)

# plot train and test-ROC curve and write down AUC values and gap between ROC-curves
plot(roc_train, col = "blue")
lines(roc_test,  col = "red")
legend("bottomright",
  legend = c(paste0("Train (AUC = ", round(auc_train_glm, 3), ")"),
    paste0("Test  (AUC = ", round(auc_test_glm, 3), ")"),
    paste0("Gap   = ", round(gap_glm, 3))),
  col = c("blue", "red", NA),lwd = c(2, 2, NA), bty = "n")

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
#---- Generalized Additive Model -----------------------------------------------
#-------------------------------------------------------------------------------

pred_cols <- intersect(filtered_var, names(goat_train))

#
var_gam <- as.formula(paste("occ ~", paste0("s(", filtered_var, ")", collapse = " + ")))

# train GAM model
model_gam <- gam(var_gam, data = goat_train, family = binomial)
summary(model_gam)

pred_train_gam <- predict(model_gam, newdata = goat_train, type = "response")
pred_test_gam  <- predict(model_gam, newdata = goat_test,  type = "response")


#---- Check for overfitting ----------------------------------------------------

roc_train_gam <- pROC::roc(goat_train$occ, pred_train_gam)
roc_test_gam  <- pROC::roc(goat_test$occ,  pred_test_gam)
auc_train_gam <- pROC::auc(roc_train_gam)
auc_test_gam  <- pROC::auc(roc_test_gam)
gap_gam       <- as.numeric(auc_train_gam - auc_test_gam)

plot(roc_train_gam, col = "blue", main = "ROC — GAM (Train vs. Test)")
lines(roc_test_gam,  col = "red")
legend("bottomright",
       legend = c(
         paste0("Train (AUC = ", round(auc_train_gam, 3), ")"),
         paste0("Test  (AUC = ", round(auc_test_gam, 3), ")"),
         paste0("Gap   = ", round(gap_gam, 3))),
       col = c("blue","red", NA), lwd = c(2,2,NA), bty = "n")


#---- Predict to raster --------------------------------------------------------

bioclim_gam <- bioclim_data[[filtered_var]]

vars_in_gam <- setdiff(all.vars(var_gam), "occ")
vars_in_gam <- intersect(vars_in_gam, names(bioclim_gam))
bioclim_gam <- bioclim_gam[[vars_in_gam]]

model_gam_pred <- terra::predict(object = bioclim_gam, model= model_gam,type = "response",na.rm = TRUE)

plot(model_gam_pred, main = "Predicted Habitat Suitability (GAM)")




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

# formular for maxnet model
fmx <- maxnet.formula(p = mx_train$occ, data = mx_train[, pred_cols], classes = "lqph")
# maxnet model 
model_maxnet <- maxnet(p= mx_train$occ, data = mx_train[, pred_cols], f= fmx)
print(model_maxnet)

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

# Random Forest model with Bootstrap-Sampling 
# saving the state of the random number generator to use an individual set.seed() in the RF model
old_seed <- .Random.seed

# same as before
pred_cols <- intersect(filtered_var, names(goat_train))

rf_train <- goat_train[, c("occ", pred_cols)]
rf_test  <- goat_test[,  c("occ", pred_cols)]

rf_train <- rf_train[complete.cases(rf_train), ]
rf_test  <- rf_test[complete.cases(rf_test), ]

rf_train$occ <- factor(rf_train$occ, levels = c(0, 1))
rf_test$occ  <- factor(rf_test$occ,  levels = c(0, 1))


set.seed(50)
model_rf <- randomForest(occ ~ ., data = rf_train, ntree = 500, importance = TRUE)
# model check:
print(model_rf)
model_rf$ntree
model_rf$mtry 
importance(model_rf)

#reset seed
.Random.seed <- old_seed


#---- Check for overfitting ----------------------------------------------------

pred_train_rf_oob <- predict(model_rf, type = "prob")[, "1"]
pred_test_rf      <- predict(model_rf, newdata = rf_test, type = "prob")[, "1"]

roc_train_rf <- pROC::roc(rf_train$occ, as.numeric(pred_train_rf_oob))
roc_test_rf  <- pROC::roc(rf_test$occ,  as.numeric(pred_test_rf))

auc_train_rf <- pROC::auc(roc_train_rf)
auc_test_rf  <- pROC::auc(roc_test_rf)
gap_rf       <- as.numeric(auc_train_rf - auc_test_rf)

plot(roc_train_rf, col = "blue", main = "Random Forest (OOB vs. Test)")
lines(roc_test_rf,  col = "red")
legend("bottomright",
       legend = c(
         paste0("Train/OOB (AUC = ", round(auc_train_rf, 3), ")"),
         paste0("Test (AUC = ", round(auc_test_rf, 3), ")"),
         paste0("Gap= ", round(gap_rf, 3))),
       col = c("blue","red", NA), lwd = c(2,2,NA), bty = "n")


#---- Predict to raster --------------------------------------------------------

bioclim_rf <- bioclim_data[[pred_cols]]

model_rf_pred <- terra::predict(bioclim_rf, model_rf, type = "prob", index = 2, na.rm = TRUE)
plot(model_rf_pred, main = "Predicted Habitat Suitability (Random Forest)")




################################################################################
#-------------------------------------------------------------------------------
#---- Ensemble Model -----------------------------------------------------------
#-------------------------------------------------------------------------------

# mean value for all predicted test and train data
ens_test <- rowMeans(cbind(
  as.numeric(test_preds),
  as.numeric(pred_test_gam),
  as.numeric(pred_test_mx),
  as.numeric(pred_test_rf)
), na.rm = TRUE)

ens_train <- rowMeans(cbind(
  as.numeric(train_preds),
  as.numeric(pred_train_gam),
  as.numeric(pred_train_mx),
  as.numeric(pred_train_rf_oob)
), na.rm = TRUE)

# Roc, Auc and Gap
roc_train_ens <- roc(goat_train$occ, ens_train)
roc_test_ens  <- roc(goat_test$occ,  ens_test)

auc_train_ens <- auc(roc_train_ens)
auc_test_ens  <- auc(roc_test_ens)
gap_ens       <- as.numeric(auc_train_ens - auc_test_ens)

plot(roc_train_ens, col="blue", main="ROC — Ensemble (Train vs. Test)")
lines(roc_test_ens, col="red")
legend("bottomright",
       legend = c(
         paste0("Train (AUC = ", round(auc_train_ens, 3), ")"),
         paste0("Test  (AUC = ", round(auc_test_ens, 3), ")"),
         paste0("Gap   = ", round(gap_ens, 3))
       ),
       col=c("blue","red",NA), lwd=c(2,2,NA), bty="n")


#---- Predict to raster --------------------------------------------------------
#mean value of all predictions
model_ens_pred <- (model_glm_pred + model_gam_pred + model_maxnet_pred + model_rf_pred) / 4

plot(model_ens_pred, main = "Predicted Habitat Suitability, Mean Ensemble")




# ################################################################################
# #-------------------------------------------------------------------------------
# #---- Compare Models -----------------------------------------------------------
# #-------------------------------------------------------------------------------
# # PROBLEM: Cate has no methods for MaxEnt
# # Prep training data
# cv_data <- goat[, c("occ", filtered_var)]
# cv_data <- cv_data[complete.cases(cv_data), ]
# 
# # caret requires labels to be factors with levels
# cv_data$occ <- factor(ifelse(cv_data$occ == 1, "Presence", "Absence"), levels = c("Presence", "Absence"))
# 
# # cross validation (10-times (collect AUC-values for 10 folds)) 
# ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = "final")
# 
# old_seed <- .Random.seed
# set.seed(50)
# 
# models_cv <- list(
#   glm = train(occ ~ ., data = cv_data, method = "glm",
#               family = binomial, trControl = ctrl, metric = "ROC"),
#   gam = train(occ ~ ., data = cv_data, method = "gam",
#               trControl = ctrl, metric = "ROC"),
#   rf  = train(occ ~ ., data = cv_data, method = "rf",
#               trControl = ctrl, metric = "ROC"))
# 
# .Random.seed <- old_seed
# 
# # saving metrics (AUC-Values) for every train fold of every model)
# resamps <- caret::resamples(models_cv)
# # compare boxplots of AUCs
# bwplot(resamps, metric = "ROC", main = "Model Comparison (AUC)")
# 
# 
################################################################################
#-------------------------------------------------------------------------------
#---- Threshold Model ----------------------------------------------------------
#-------------------------------------------------------------------------------
# Representation of the habitat when all grid cells with a certain threshold are
# combined from all models. So, the locations where all models predict a very
# high probability of occurrence of the Mountain Goat.
# The condition for identifying a grid cell as a potential habitat is that at
# least three of the four models have a prediction value of >= 0.75.

# threshold
thr <- 0.70

# stacking the four main models
preds <- c(model_glm_pred, model_gam_pred, model_rf_pred, model_maxnet_pred)
names(preds) <- c("GLM","GAM","RF","MaxEnt")

# count cells with valid value
hits <- app(preds, fun = function(x) sum(x >= thr, na.rm = TRUE))
habitat_bin <- ifel(hits >= 3, 1, 0)

levels(habitat_bin) <- data.frame(ID = c(0,1), class = c("Probably no habitat", "Habitat"))

plot(habitat_bin, main = "Consensus Habitat (≥3 models with predicted value ≥0.75)",
     col = c("grey80", "darkgreen"), legend = FALSE)
legend("bottomleft",
       legend = c("Probably no habitat", "Habitat"),
       fill   = c("grey80", "darkgreen"),
       bty    = "n")

#### Kicking out the GLM model

# stacking the four main models
preds_without_glm <- c( model_gam_pred, model_rf_pred, model_maxnet_pred)
names(preds_without_glm) <- c("GAM","RF","MaxEnt")

# count cells with valid value
hits_2 <- app(preds_without_glm, fun = function(x) sum(x >= thr, na.rm = TRUE))
habitat_bin_2 <- ifel(hits >= 2, 1, 0)

levels(habitat_bin_2) <- data.frame(ID = c(0,1), class = c("Probably no habitat", "Habitat"))

plot(habitat_bin_2, main = "Consensus Habitat (≥2 models with predicted value ≥0.75)",
     col = c("grey80", "darkgreen"), legend = FALSE)
plot(st_geometry(world_sf), add = TRUE, col = "transparent", border = "grey50", lwd = 0.1)
legend("bottomleft",
       legend = c("Probably no habitat", "Habitat"),
       fill   = c("grey80", "darkgreen"),
       bty    = "n")




################################################################################
##         All predicted rasters                                              ##
################################################################################

plot_suffix <- if (elev) "with_elev" else "without_elev"
plot_file   <- paste0("Suitability_models_", plot_suffix, ".png")

png(filename = plot_file, width = 2000, height = 1500, res = 150)

par(mfrow = c(2, 3))
plot(model_glm_pred, main = "GLM Habitat Suitability")
plot(model_gam_pred, main = "GAM Habitat Suitability")
plot(model_maxnet_pred, main = "MaxEnt Habitat Suitability")
plot(model_rf_pred, main = "Random Forest Suitability")
plot(model_ens_pred, main = "Predicted Habitat Suitability, Mean Ensemble")
par(mfrow = c(1, 1))

dev.off()
cat("Plot got saved as:", plot_file, "\n")

################################################################################
##         Table with all AUC Values                                          ##
################################################################################

results <- data.frame(
  Model     = c("GLM", "GAM", "MaxEnt", "Random Forest", "Ensemble"),
  AUC_train = round(c(
    as.numeric(auc_train_glm),
    as.numeric(auc_train_gam),
    as.numeric(auc_train_mx),
    as.numeric(auc_train_rf),
    as.numeric(auc_train_ens)
  ), 4),
  AUC_test  = round(c(
    as.numeric(auc_test_glm),
    as.numeric(auc_test_gam),
    as.numeric(auc_test_mx),
    as.numeric(auc_test_rf),
    as.numeric(auc_test_ens)
  ), 4)
)

results$GAP <- round(results$AUC_train - results$AUC_test, 4)

# creating table depending if elev = TRUE or FALSE
if (elev) {
  results_with_elev <- results
} else {
  results_without_elev <- results}

# individual naming
file_suffix <- if (elev) "with_elev" else "without_elev"
file_name   <- paste0("AUC_results_", file_suffix, ".xlsx")

# create .xlsx table in working directory
write.xlsx(results, file = file_name, rowNames = FALSE)
cat("Table got saved as:", file_name, "\n")


# TODO: GLCM
# TODO: Plot achsenbeschriftung long und lat, legendenbeschriftung, irgendwie richtigen plot machen
# TODO: besseren threshold plot auch noch hinzufügen
# TODO: interpret model summary
