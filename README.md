# Summer Term 2025 | Spatial Modelling and Prediction | 04-GEO-OMA1

This habitat modelling was created as part of the *Spatial Modelling and Prediction course* held at the Earth Observation Research Cluster.

R version: 4.4.1

# Data
**Species Data:** Global Core Biodata Resource (2025): *Oreamnos americanus (Blainville, 1816)*, URL: https://doi.org/10.15468/39omei

**Bioclimatic & Eleveation Variables:** Fick, S. E., & Hijmans, R. J. (2017): *WorldClim 2: new 1km spatial resolution climate surfaces for global land areas.*, URL: https://www.worldclim.org/data/worldclim21.html

# Content

The script models the potential habitat of the **_Oreamnos americanus_** using a **Generalized Liear Model, Generalized Additive Model, Maximum Entropy Model**, and **Random Forest** model. The **_Oreamnos americanus_** distribution is highly fragmented into niches, it has been introduced to some areas and naturally inhabits mountain regions.

### Packages
```r
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
```

### 1. Data download of the species
To download the complete species data, you need your own account (not a Google account). Your username, password and email address must be entered in the download code. 
```r
path <- setwd("...")

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
```

### 2. Filtering the data so that only animals living in freedom are included in the model.
```r
# Keeping only entries with accuracy of 1 km, and deleting entries of preserved specimen
occ_data_filtered <- occ_data %>%
  filter(
    (is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= 1000) &
      !basisOfRecord %in% c("PRESERVED_SPECIMEN")&
      decimalLongitude < -62)

# remove unnecessary columns
occ_data_filtered <- occ_data_filtered %>% select(c("decimalLongitude","decimalLatitude"))
```

### 3. Defining AOI based on observations with distance from the edge.
```r
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
```
```r
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
```
<img width="565" height="368" alt="AOI" src="https://github.com/user-attachments/assets/011862c1-cd99-495d-b296-60ca72a58681" />

### 4. Data download of bioclimatic and elevation variables and cropping to the AOI

```r
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
```
### 5. Creating 5000 pseudo points in AOI (except ocean) and linking them to the species dataset. The bioclimatic and altitude variables are assigned to the points and only one point is retained per grid cell.

```r
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
```
```r
# plot observation and random seed points
ggplot() +
  geom_sf(data = world_sf, color = "grey50", fill = "white", linewidth = 0.2) +
  geom_point(data = goat, aes(x = decimalLongitude, y = decimalLatitude, color = factor(occ)),size = 0.5) +
  scale_color_manual(name= "Legend", values = c("0" = "grey30", "1" = "red"), labels = c("Random", "Observed Occurrence")) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  labs(title = "Observed occurrence & Random",
       x= "Longitude",
       y= "Latitude")
```
<img width="572" height="368" alt="AOI_random" src="https://github.com/user-attachments/assets/5be057a6-b3ee-426d-888b-0e18aca97b1c" />

### 6. After trying out various combinations of variables, highly correlated variables ( cor. >=0.7 / -0.7) were removed from the data set. 
```r
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
```
The elevation variable can be used separately: At this point in the code, the elevation variable can be included in the model calculation by setting *elev* to *TRUE*.
```r
  elev <- FALSE 

  if (elev) {
    filtered_var <- vars_with_elev
  } else {
    filtered_var <- vars_without_elev
  }
```

### 7. The AOI was divided into hexagons with a width of 200 km, and the data was divided into train or test according to the hexagon in order to reduce spatial autocorrelation. This train-test division was randomly assigned 10 times (10-folds).

```r
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
  plot = TRUE)

.Random.seed <- old_seed

folds <- cv$folds_list

```
<img width="1260" height="1008" alt="hex" src="https://github.com/user-attachments/assets/89133a94-57f3-4057-ad6f-d6d7f5f52f9d" />

### 8. The four models are trained in a large loop with each of the 10 folds, and their mean prediction is calculated and stored. The mean AUC and AUC standard deviation are also determined in this way. In addition, an ensemble model is generated that determines the mean value from all model predictions. 

Preparation for storing intermediate results of the folds.
```r
# result dataframe for all folds
results_all <- data.frame()

# create list to save prediction for every fold prediction of every model
mean_preds <- list(GLM=NULL, GAM=NULL, MaxNet=NULL, RF=NULL, Ensemble=NULL)

r_all <- bioclim_data[[filtered_var]]
```

!!! The loop may take several minutes to complete !!!
```r
for (i in seq_along(folds)) {
  # outputs the current fold to show progress
  cat("\n--- Fold", i, "---\n")
  
#---- create train and test data -----------------------------------------------
  train_idx <- unlist(folds[[i]][[1]])
  test_idx  <- unlist(folds[[i]][[2]])
  
  goat_train <- goat[train_idx, ]
  goat_test  <- goat[test_idx, ]
```

#### Generalized Linear Model
```r
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
```

#### Generalized Additive Model
```r
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
```

#### Maximum Entropy Model
```r
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

#---- Predict to raster --------------------------------------------------------  
  # extract relevant predictors
  r_mx <- r_all[[pred_cols]]
  
  # MaxEnt model using for prediction for each raster zell (Values between 0 and 1 (cloglog))
  pred_r_mx <- terra::predict(r_mx, model_maxnet, type = "cloglog", na.rm = TRUE)
```

#### Random Forest Model
```r
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
```

#### Ensemble Model
```r
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
```

#### Saving all results
```r
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
```
### 9. Plots and tables were saved in the working directory.

#### Table with all AUC Values
```r
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
```
```r
# creating table depending if elev = TRUE or FALSE
file_suffix <- if (elev) "with_elev" else "without_elev"
# individual naming
file_name   <- paste0("AUC_results_", file_suffix, ".xlsx")

# create .xlsx table in working directory
openxlsx::write.xlsx(summary_results_rounded, file = file_name, rowNames = FALSE)
cat("Table got saved as:", file_name, "\n")
```

#### All predicted rasters

```r
# calculating the average prediction raster for each model
k_folds <- length(folds)

for (m in names(mean_preds)) {
  mean_preds[[m]] <- mean_preds[[m]] / k_folds}
```

Plot
```r
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
```

### 10. A threshold model with a binary grid was created and shows which cells from at least three models have a probability of at least 70% of being the habitat of the species.
```r
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

```
### 11. A standard deviation grid and histogram of all model results and a spatial correlation matrix were created.
#### Grid of the standard deviation of the predictions of all models 
```r
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
```

#### Histogram of the standard deviation
```r
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
```

Correlation matrix of the spatial prediction of the models for the potential habitat, using a sample size of 50,000 cells 
```r
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
```


# Example results

## AUC Table

AUC results without elevation variable
| Model    | Mean_AUC_train | SD_AUC_train | Mean_AUC_test | SD_AUC_test | Mean_Gap | SD_Gap |
|----------|---------------:|-------------:|--------------:|------------:|---------:|-------:|
| Ensemble | 0.9188         | 0.0068       | 0.8969        | 0.0643      | 0.0219   | 0.0663 |
| GAM      | 0.9273         | 0.0064       | 0.8468        | 0.0730      | 0.0805   | 0.0778 |
| GLM      | 0.8691         | 0.0085       | 0.8592        | 0.0792      | 0.0099   | 0.0864 |
| MaxNet   | 0.9257         | 0.0081       | 0.8710        | 0.0825      | 0.0547   | 0.0891 |
| RF       | 0.9582         | 0.0050       | 0.8975        | 0.0621      | 0.0607   | 0.0655 |


AUC results with elevation variable
| Model    | Mean_AUC_train | SD_AUC_train | Mean_AUC_test | SD_AUC_test | Mean_Gap | SD_Gap |
|----------|---------------:|-------------:|--------------:|------------:|---------:|-------:|
| Ensemble | 0.9285         | 0.0063       | 0.9143        | 0.0484      | 0.0142   | 0.0486 |
| GAM      | 0.9352         | 0.0056       | 0.8606        | 0.0634      | 0.0746   | 0.0677 |
| GLM      | 0.8697         | 0.0083       | 0.8522        | 0.0835      | 0.0175   | 0.0904 |
| MaxNet   | 0.9359         | 0.0068       | 0.8812        | 0.0676      | 0.0548   | 0.0729 |
| RF       | 0.9715         | 0.0031       | 0.9275        | 0.0482      | 0.0440   | 0.0501 |

The model with elevation data generally has a slightly higher AUC value, and all models show good discriminative power, ranking a presence point higher than an absence point. The low values for Mean_Gap are a small overfitting test. Values between 0.05 and 0.01 are good and indicate that the model works a little better with the training data, but that there is probably no overfitting. Values above 0.05, as is the case with GAM and MaxNet, indicate a tendency towards overfitting.

## Mean Model Prediction

Similar patterns appear to exist, but the predictions of MaxNet and GAM (both of which may tend to overfit) correlate more strongly with each other, as well as MaxNet and GLM. The spatial correlation matrices are also shown at the end of the readme file.
Prediction grids were created based on the average prediction (10-fold runs) of the models.

**GLM:** is a linear model, which is why it tends to produce smooth gradients and less capture of niches --> (See maps: large areas with values in the middle range).

**GAM:** Non-linear model and shows local hotspots in the maps. Possibly reduce degrees of freedom, as there is a tendency to overfit.

**MaxNet:** Generally higher values and strongly separated, large niches.

**RF:** Can handle strong fragmentation well by detecting small-scale patterns, which is also reflected in the maps.


Without elevation variable

<img width="3600" height="1500" alt="Suitability_models_mean_without_elev" src="https://github.com/user-attachments/assets/d81b9ca5-abc4-4751-933b-3b30e475c259" />

With elevation variable. 

<img width="3600" height="1500" alt="Suitability_models_mean_with_elev" src="https://github.com/user-attachments/assets/df8139e0-aec9-4bb1-a477-ab85183ee3dd" />

## Models Threshold

The fact that the models tend to disagree is shown by the low overlap of the high habitat values.
Consistent model prediction for habitat suitability probability values => 0.7.
Without elevation variable.

<img width="2000" height="1500" alt="Threshold_models_mean_without_elev" src="https://github.com/user-attachments/assets/5afb3c0a-7fdf-4923-be64-3b3d52578a1c" />

With elevation variable

<img width="2000" height="1500" alt="Threshold_models_mean_with_elev" src="https://github.com/user-attachments/assets/83afd1a0-8bc2-42dc-951c-6aa5d242a2de" />

## Standard deviation map and histogram of model predictions

The standard deviation measure also reflects the uncertainties in the more precise habitat determination. 0 = complete agreement, 0.5 = high prediction conflicts (two models say yes, two say no).

Without elevation variable.
<p align="center">
  <img src="https://github.com/user-attachments/assets/1d565232-e8a3-42aa-ba41-f78420c3f245"
       alt="SD_between_models_without_elev"
       width="48%" />
  <img src="https://github.com/user-attachments/assets/c1b04bd2-3641-4117-9fb0-c2b47487771f"
       alt="SD_histogram_without_elev"
       width="48%" />
</p>

With elevation variable

<p align="center">
  <img src="https://github.com/user-attachments/assets/55bc0108-9f28-4343-ae67-c6ae08539e24"
       alt="SD_between_models_with_elev"
       width="48%" />
  <img src="https://github.com/user-attachments/assets/59af131b-393e-4c4e-b40d-f4959452cf21"
       alt="SD_histogram_with_elev"
       width="48%" />
</p>


## Model prediction correlation 

Spatial correlation of the models, how similar are their predictions of the probability of the habitat occurring?

<p align="center">
  <img src="https://github.com/user-attachments/assets/5346b887-6327-47c2-833d-9d5a711153ed"
       alt="Model_prediction_correlation_without_elev"
       width="48%" />
  <img src="https://github.com/user-attachments/assets/f12ec48b-881a-4463-a4d2-030a2a3c1cdc"
       alt="Model_prediction_correlation_with_elev"
       width="48%" />
</p>
<p align="center">
  Without elevation variable &nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp; With elevation variable
</p>

The models that include height do not necessarily have greater overlap. The performance of the models depends on the randomly distributed train test data set.
# Challenges

_The species has a very specific, fragmented habitat in mountainous areas. With a finer grid (currently approx. 4 km resolution), the differences could be better captured, but this requires a lot of computing power. 

_Initially, the model was trained with a classic division of points in a 70/30 ratio, which led to better AUC values, but the blockCV package was then used to prevent spatial autocorrelation.

_Due to the small habitats and the condition that there can only be one point per cell, over 3,000 of the 4,177 presence observations were deleted after filtering (occ_data_filtered). This resulted in a background to presence ratio of 1:7.

_The inclusion of strongly correlated bioclimatic variables in various combinations resulted in AUC values between 0.95 and 0.99. After filtering these out, more specific variables remained, which are also used in the current script version (Bio 8, 10, 13, 15, 18,(elevation)).
