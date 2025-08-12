################################################################################
#### Habitat modeling - Oreamnos americanus ####################################
################################################################################

library(rgbif)
library(dplyr)
library(maps)
library(terra)
library(geodata)
library(sf)
library(ggplot2)

#---- Download data of "Oreamnos americanus (Blainville, 1816)" ----------------

path <- setwd("C:/users/Duck/Documents/Studium/EAGLE/2_semester/3_Spatio_Temporal_Modelling/ST25_SpatioTemporalModeling_Exam/ST25_SpatioTemporalModeling/data")

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

# TODO: Need to fic path
occ_data <- occ_download_import("0003051-250811113504898.zip")


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
xmin <- -170   # West Alasca
xmax <-  -52   # East Canada
ymin <-   33   # Somewhere South USA
ymax <-   68   # Somewhere North canada

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

#---- ----
