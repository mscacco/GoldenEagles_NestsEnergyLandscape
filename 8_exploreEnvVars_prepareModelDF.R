#______________________________________________________
# Explore ENV predictors and prepare final Dataset ####
#______________________________________________________

# This script the difference in environmental variables obtained at different scales and resolutions
# Explore their distribution relative to the nest presences and absences
# Exports a clean dataset with relevant variables ready for the models

library(sf)
library(ggplot2)

setwd("/home/martina/ownCloud/Martina/ProgettiVari/GoldenEagles_COSMO/Nests_EnergyLandscape_Giulia/Giulia_sharedFolder/")

#this need to be updated once i have downloaded tiles
model_data <- readRDS("outputData/PCA_wind_pred_2708cio_allCRS_lulc_topo.rds")


#____________________
# Explore env variables at different resolutions


#____________________
# Explore missing values in DEM

# They are only absenses and they are at the margin of the swiss territory.
library(mapview)
mapview(dem_mean_100) +
  mapview(model_data[is.na(model_data$dem_2m), ], col.regions = "red")

table(model_data$response[is.na(model_data$dem_2m)])
terr <- table(model_data$territory_ID[is.na(model_data$dem_2m)])
(terr/table(model_data$territory_ID[model_data$territory_ID %in% names(terr)]))*100

# for 4 territories, the number of absences with NA values is > 20% of the N of locations
# so either we remove those 4 territories, or we use a 30m aster dem, or we replace only the missing values with the aster dem

4/length(unique(model_data$territory_ID))
plot(model_data$dtm, model_data$dem_2m)
plot(model_data$dtm, model_data$dem_mean_100m)
plot(model_data$slope, model_data$slope_2m)
plot(model_data$slope, model_data$slope_100m)


#___________________
# Compare variables between nests presence/absence

# Topography
ggplot(model_data, aes(x = dem_2m, colour = factor(response))) +
  geom_density(linewidth = 1) +
  labs(x = "Terrain elevation (2 m pixels)",y = "Density",colour = "Response") +
  theme_minimal()
ggplot(model_data, aes(x = slope_max_100m, colour = factor(response))) +
  geom_density(linewidth = 1) +
  labs(x = "Maximum terrain slope (100 m pixels)",y = "Density",colour = "Response") +
  theme_minimal()
ggplot(model_data, aes(x = slope_100m, colour = factor(response))) +
  geom_density(linewidth = 1) +
  labs(x = "Terrain slope (100 m pixels)",y = "Density",colour = "Response") +
  theme_minimal()
ggplot(model_data, aes(x = aspect_100m, colour = factor(response))) +
  geom_density(linewidth = 1) +
  labs(x = "Terrain aspect (100 m pixels)",y = "Density",colour = "Response") +
  theme_minimal()

# Orographic Uplift
ggplot(model_data, aes(x = OroUplift_100m, colour = factor(response))) +
  geom_density(linewidth = 1) +
  labs(x = "Orographic uplift",y = "Density",colour = "Response") +
  theme_minimal()

# Vertical distance
ggplot(model_data, aes(x = verticalDist_nest_fromMedTerrElev, colour = factor(response))) +
  geom_density(linewidth = 1) +
  labs(x = "Nest - median territory elevation",y = "Density",colour = "Response") +
  theme_minimal()




