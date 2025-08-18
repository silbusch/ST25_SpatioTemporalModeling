# Summer Term 2025 | Spatial Modelling and Prediction | 04-GEO-OMA1

This habitat modelling was created as part of the *Spatial Modelling and Prediction course* held at the Earth Observation Research Cluster.

R version: 4.4.1

## Data
**Species Data:** Global Core Biodata Resource (2025): *Oreamnos americanus (Blainville, 1816)*, URL: https://doi.org/10.15468/39omei

**Bioclimatic & Eleveation Variables:** Fick, S. E., & Hijmans, R. J. (2017): *WorldClim 2: new 1km spatial resolution climate surfaces for global land areas.*, URL: https://www.worldclim.org/data/worldclim21.html

## Content

The script models the potential habitat of the **_Oreamnos americanus_** using a **Generalized Liear Model, Generalized Additive Model, Maximum Entropy Model**, and **Random Forest** model. The **_Oreamnos americanus_** distribution is highly fragmented into niches, it has been introduced to some areas and naturally inhabits mountain regions.

#### 1. Data download of the species
To download the complete species data, you need your own account (not a Google account). Your username, password and email address must be entered in the download code. The following is a brief description of the procedure and some code snippets for preparation analysis are shown. The script itself contains further comments on the individual steps.

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
<img width="565" height="368" alt="AOI" src="https://github.com/user-attachments/assets/011862c1-cd99-495d-b296-60ca72a58681" />

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

<img width="572" height="368" alt="AOI_random" src="https://github.com/user-attachments/assets/5be057a6-b3ee-426d-888b-0e18aca97b1c" />

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

```r
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

folds <- cv$folds_list
```
<img width="1260" height="1008" alt="hex" src="https://github.com/user-attachments/assets/89133a94-57f3-4057-ad6f-d6d7f5f52f9d" />

#### 8. The four models are trained in a large loop with each of the 10 folds, and their mean prediction is calculated and stored. The mean AUC and AUC standard deviation are also determined in this way. In addition, an ensemble model is generated that determines the mean value from all model predictions.
#### 9. Plots and tables were saved in the working directory.
#### 10. A threshold model with a binary grid was created and shows which cells from at least three models have a probability of at least 70% of being the habitat of the species.
#### 11. A standard deviation grid of all model results and a spatial correlation matrix were created.



## Example results

### AUC Table

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

### Mean Model Prediction

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

### Models Threshold

The fact that the models tend to disagree is shown by the low overlap of the high habitat values.
Consistent model prediction for habitat suitability probability values => 0.7.
Without elevation variable.

<img width="2000" height="1500" alt="Threshold_models_mean_without_elev" src="https://github.com/user-attachments/assets/5afb3c0a-7fdf-4923-be64-3b3d52578a1c" />

With elevation variable

<img width="2000" height="1500" alt="Threshold_models_mean_with_elev" src="https://github.com/user-attachments/assets/83afd1a0-8bc2-42dc-951c-6aa5d242a2de" />

### Standard deviation map and histogram of model predictions

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


### Model prediction correlation 

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
## Challenges

_The species has a very specific, fragmented habitat in mountainous areas. With a finer grid (currently approx. 4 km resolution), the differences could be better captured, but this requires a lot of computing power. 

_Initially, the model was trained with a classic division of points in a 70/30 ratio, which led to better AUC values, but the blockCV package was then used to prevent spatial autocorrelation.

_Due to the small habitats and the condition that there can only be one point per cell, over 3,000 of the 4,177 presence observations were deleted after filtering (occ_data_filtered). This resulted in a background to presence ratio of 1:7.

_The inclusion of strongly correlated bioclimatic variables in various combinations resulted in AUC values between 0.95 and 0.99. After filtering these out, more specific variables remained, which are also used in the current script version (Bio 8, 10, 13, 15, 18,(elevation)).
