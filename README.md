# Summer Term 2025 | Spatial Modelling and Prediction | 04-GEO-OMA1

This habitat modelling was created as part of the *Spatial Modelling and Prediction course* held at the Earth Observation Research Cluster.

R version: 4.4.1


## Content

The script models the potential habitat of the **_Oreamnos americanus_** using a Generalized Liear Model, Generalized Additive Model, Maximum Entropy Model, and Random Forest model. 

#### 1. Data download of the species
To download the complete species data, you need your own account (not a Google account). Your username, password and email address must be entered in the download code. The following is a brief description of the procedure and some corde snippets are shown. The script itself contains further comments on the individual steps.

```r
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

#### 2. Filtering the data so that only animals living in freedom are included in the model (e.g. no preserved specimen).
#### 3. Defining AOI based on observations with buffer
<img width="447" height="368" alt="5e067d0b-8242-468b-9bd7-1e6e36639d7d" src="https://github.com/user-attachments/assets/c63358e8-f455-46b1-9647-a66ce21d0c44" />

#### 4. Data download of bioclimatic and elevation variables and cropping to the AOI

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
#### 5. Creating 5000 pseudo points in AOI (except ocean) and linking them to the species dataset. The bioclimatic and altitude variables are assigned to the points and only one point is retained per grid cell.

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
<img width="523" height="368" alt="976d364d-73b8-4018-98b7-70bdcb348e79" src="https://github.com/user-attachments/assets/de8c9a40-0d20-472a-bd8d-1fa12927ee94" />

#### 6. After trying out various combinations of variables, highly correlated variables ( cor. >=0.7 / -0.7) were removed from the data set. The elevation variable can be used separately.

At this point in the code, the elevation variable can be included in the model calculation by setting *elev* to *TRUE*.
```r
  elev <- FALSE 

  if (elev) {
    filtered_var <- vars_with_elev
  } else {
    filtered_var <- vars_without_elev
  }
```

#### 7. The AOI was divided into hexagons with a width of 200 km, and the data was divided into train or test according to the hexagon in order to reduce spatial autocorrelation. This train-test division was randomly assigned 10 times (10-folds).

<img width="1920" height="1008" alt="hex" src="https://github.com/user-attachments/assets/cbb590be-6c8e-4d9a-87ec-36e1ac85757a" />

#### 8. The four models are trained in a large loop with each of the 10 folds, and their mean prediction is calculated and stored. The mean AUC and AUC standard deviation are also determined in this way. In addition, an ensemble model is generated that determines the mean value from all model predictions.
#### 9. Plots and tables were saved in the working directory.
#### 10. A threshold model with a binary grid was created and shows which cells from at least three models have a probability of at least 70% of being the habitat of the species.
#### 11. A standard deviation grid of all model results and a spatial correlation matrix were created.
...
## Data


## Example results

### AUC Table

AUC results without elevation variable
| Model       | Mean_AUC_train | SD_AUC_train  | Mean_AUC_test | SD_AUC_test    | Mean_Gap      | SD_Gap  |
| :-----------| :-----------: | -------------: |-------------: |--------------: |------------: |---------: |
| Ensemble    | 0,9188        | 0,0068         | 0,8969        | 0,0643         | 0,0219       | 0,0663
| GAM         | 0,9273        | 0,0064         | 0,8468        | 0,073          | 0,0805       | 0,0778
| GLM         | 0,8691        | 0,0085         | 0,8592        | 0,0792         | 0,0099       | 0,0864
| MaxNet      | 0,9257        | 0,0081         | 0,871         | 0,0825         | 0,0547       | 0,0891
| RF          | 0,9582        | 0,005          | 0,8975        | 0,0621         | 0,0607       | 0,0655

AUC results with elevation variable
| Model       | Mean_AUC_train | SD_AUC_train  | Mean_AUC_test | SD_AUC_test    | Mean_Gap      | SD_Gap  |
| :-----------| :-----------: | -------------: |-------------: |--------------: |------------: |---------: |
| Ensemble    | 0,9285        | 0,0063         | 0,9143        | 0,0484         | 0,0142       | 0,0486
| GAM         | 0,9352        | 0,0056         | 0,8606        | 0,0634          | 0,0746       | 0,0677
| GLM         | 0,8697        | 0,0083         | 0,8522        | 0,0835         | 0,0175       | 0,0904
| MaxNet      | 0,9359        | 0,0068         | 0,8812         | 0,0676         | 0,0548       | 0,0729
| RF          | 0,9715        | 0,0031          | 0,9275        | 0,0482         | 0,044       | 0,0501

### Mean Model Prediction

Without elevation variable

<img width="2000" height="1500" alt="Suitability_models_mean_without_elev" src="https://github.com/user-attachments/assets/4fece05b-cdad-46bf-a394-aa3f8c9d7180" />

With elevation variable

<img width="2000" height="1500" alt="Suitability_models_mean_with_elev" src="https://github.com/user-attachments/assets/31066a68-5f57-436c-bb55-e55215fdf511" />

### Models Threshold
Consistent model prediction for habitat suitability probability values => 0.7.
Without elevation variable.

<img width="2000" height="1500" alt="Threshold_models_mean_without_elev" src="https://github.com/user-attachments/assets/c9891d02-3798-4711-86d8-c4a2af0faa0d" />

With elevation variable

<img width="2000" height="1500" alt="Threshold_models_mean_with_elev" src="https://github.com/user-attachments/assets/21bec3e9-49c4-48b7-ac01-098428804b57" />

### Standard deviation map and histogram of model predictions

The standard deviation between the model predictions. 0 = complete agreement, 0.5 = high prediction conflicts (two models say yes, two say no).

Without elevation variable.

<p align="center">
  <img src="https://github.com/user-attachments/assets/750dbbe9-36b5-4ed0-8f48-a3ea86f16630"
       alt="SD_between_models_without_elev"
       width="48%" />
  <img src="https://github.com/user-attachments/assets/626aeb20-9ce4-4176-8935-031a9ca82cdb"
       alt="SD_histogram_without_elev"
       width="48%" />
</p>

With elevation variable

<p align="center">
  <img src="https://github.com/user-attachments/assets/96b125c1-aa8b-40a4-8595-d99234158ac0"
       alt="SD_between_models_with_elev"
       width="48%" />
  <img src="https://github.com/user-attachments/assets/1d292b02-b837-47dc-affe-efe34a0d413b"
       alt="SD_histogram_with_elev"
       width="48%" />
</p>


### Model prediction correlation 

Spatial correlation of the models, how similar are their predictions of the probability of the habitat occurring?

<p align="center">
  <img src="https://github.com/user-attachments/assets/60a5528a-93d7-452a-9dca-371c3e827927"
       alt="Model_prediction_correlation_without_elev"
       width="48%" />
  <img src="https://github.com/user-attachments/assets/579b568d-4075-4875-baae-d17c9c5823c6"
       alt="Model_prediction_correlation_with_elev"
       width="48%" />
</p>
<p align="center">
  Without elevation variable &nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp; With elevation variable
</p>



