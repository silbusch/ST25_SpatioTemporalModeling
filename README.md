# Summer Term 2025 | Spatial Modelling and Prediction | 04-GEO-OMA1

.....

## Data

...

## Todo on the first script run

### Data Download

### Variable decision

```r
  elev <- FALSE 

  if (elev) {
    filtered_var <- vars_with_elev
  } else {
    filtered_var <- vars_without_elev
  }
```


AUC results without elevation variable

| Model       | Mean_AUC_train | SD_AUC_train  | Mean_AUC_test | SD_AUC_test    | Mean_Gap      | SD_Gap  |
| :-----------| :-----------: | -------------: |-------------: |--------------: |------------: |---------: |
| Ensemble    | 0,9188        | 0,0068         | 0,8969        | 0,0643         | 0,0219       | 0,0663
| GAM         | 0,9273        | 0,0064         | 0,8468        | 0,073          | 0,0805       | 0,0778
| GLM         | 0,8691        | 0,0085         | 0,8592        | 0,0792         | 0,0099       | 0,0864
| MaxNet      | 0,9257        | 0,0081         | 0,871         | 0,0825         | 0,0547       | 0,0891
| RF          | 0,9582        | 0,005          | 0,8975        | 0,0621         | 0,0607       | 0,0655
