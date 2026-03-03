#___________________________________________________
# Annotate topographic and land cover variables ####
#___________________________________________________

# This script derives and aggregates topographic layers from 2 m pixels to 200 m and extract summary information to highlight the presence of cliffs and surface structures
# It annotates both topographic and land cover data to the dataset
# computes derived variables such as orographic uplift and vertical position of nest in territory

library(sf)
library(terra)
library(ggplot2)
library(mapview)
library(dplyr)

#_____________
# Define own function to calculate orographic uplift
oroUplift <- function(u, v, slope_rad, aspect_rad) {
  # u: eastward wind component (m/s) at the surface (10 m)
  # v: northward wind component (m/s) at the surface (10 m)
  # slope_rad: angle of the slope in radians (theta)
  # aspect_rad: aspect of the slope in radians between 0 and 2pi (alpha)
  ws <- sqrt(u^2 + v^2) # wind speed at 10 m height, in m/s
  beta <- atan2(u, v) + pi # wind dir at 10 m, blowing FROM, in rad between 0 and 2pi
  Ca <- sin(slope_rad) * cos(aspect_rad - beta) #lifting coefficient (if the slope value is 0, the whole Ca coefficient will be 0, so flat slopes have no oro uplift): Ca <- sin(theta) * cos(alpha - beta)
  # Calculate the orographic uplift according to Brandes et al (2004) and Bohrer et al (2012)
  oroUpl <- ws * Ca
  return(oroUpl)
}
windDir_from_deg <- function(u, v) { 
  WD.deg <- atan2(u, v) * 180/pi  #from rad to deg as +/-180
  return((WD.deg + 180) %% 360) #reverse from blowing "to" to blowing FROM, as 0-360. The %% 360 makes sure the results stays within 0-360
}
#_____________


setwd("/home/martina/ownCloud/Martina/ProgettiVari/GoldenEagles_COSMO/Nests_EnergyLandscape_Giulia/Giulia_sharedFolder/")

# Load dataset (sf object)
model_data <- readRDS("outputData/PCA_wind_pred_2708cio_allCRS.rds")


#____________________
# Annotate land cover

lucl_10m <- rast("inputData/environmental_data/LandCover_10m_PlanComp/lulc_v02_2023.tif")
plot(lucl_10m)
(leg <- read.table("inputData/environmental_data/LandCover_10m_PlanComp/lulc_v02_legend.txt", sep=",", header = T))

# set clouds and no data to NA
lucl_10m[lucl_10m %in% c(0, 10)] <- NA
names(lucl_10m) <- "lulc_10m"
plot(lucl_10m)

# compute dominant land cover in 100 m pixels
lulc_100m <- aggregate(lucl_10m, fact = 10, fun = modal, na.rm = TRUE)
plot(lulc_100m)

# compute %trees per pixel
pct_trees_100m <- aggregate(lucl_10m, fact = 10,
  fun = function(x) mean(x == 2, na.rm = TRUE))
plot(pct_trees_100m)
# compute %built area per pixel
pct_build_100m <- aggregate(lucl_10m, fact = 10,
                            fun = function(x) mean(x == 7, na.rm = TRUE))
# compute %bare ground per pixel
pct_bare_100m <- aggregate(lucl_10m, fact = 10,
                            fun = function(x) mean(x == 8, na.rm = TRUE))

# save new layers
lulc_stack_100m <- c(lulc_100m, pct_trees_100m, pct_build_100m, pct_bare_100m)
names(lulc_stack_100m) <- c("lulc_100m", "pct_trees_100m", "pct_build_100m", "pct_bare_100m")
writeRaster(lulc_stack_100m, "inputData/environmental_data/LandCover_10m_PlanComp/lulc_derived100m_pcts.tif")

# annotate data with info from all layers
model_data_prj <- st_transform(model_data, crs(lucl_10m))
val_lulc_10m <- terra::extract(lucl_10m, model_data_prj, method="nearest")
val_lulc_100m <- terra::extract(lulc_stack_100m, model_data_prj, method="nearest")

model_data_lulc <- cbind(model_data, lulc_10m=val_lulc_10m[,-1], val_lulc_100m[,-1])
table(model_data_lulc$lulc_10m)
table(model_data_lulc$lulc_100m) # they are very similar, we retain the one at 100m


# add legend and save
length(unique(model_data_lulc$ID_Tom))==nrow(model_data_lulc)
model_data_lulc <- merge(model_data_lulc, leg, by.x="lulc_100m", by.y="value")
model_data_lulc <- dplyr::rename(model_data_lulc, lcClass_100m = class)
table(model_data_lulc$lulc_100m)
table(model_data_lulc$lcClass_100m)

saveRDS(model_data_lulc, "outputData/PCA_wind_pred_2708cio_allCRS_lulc.rds")


#____________________________________
# Derive and annotate topography ####

# import the data annotated with the land cover
model_data <- readRDS("outputData/PCA_wind_pred_2708cio_allCRS_lulc.rds")

# load all tiles from HD to R as a virtual raster (does NOT load all into memory)
output_fld <- "/media/martina/Intenso/DEM_Swisstopo_Grisons_2m/"
tiles <- list.files(output_fld, pattern = "\\.tif", full.names = TRUE)
dem_vrt <- terra::vrt(tiles)

# Calculate 4 complementary topographic metrics: 
# Broad-scale slope (1) and aspect (2), calculated from a 100 m DEM obtained by aggregating the original 2 m DEM using the mean
dem_mean_100 <- terra::aggregate(dem_vrt, fact = 50, fun = mean, na.rm = TRUE)
slope_100 <- terra::terrain(dem_mean_100, v = "slope", unit = "degrees")
aspect_100 <- terra::terrain(dem_mean_100, v = "aspect", unit = "degrees")

# Fine-scale terrain ruggdeness, calculated as the maximum dem (3) and maximum slope (4) within 100 m windows based on dem at 2 m and slope computed at 2 m resolution.
slope_2m <- terra::terrain(dem_vrt, v = "slope", unit = "degrees")
slope_max_100 <- terra::aggregate(slope_2m, fact = 50, fun = max, na.rm = TRUE)
dem_max_100 <- terra::aggregate(dem_vrt, fact = 50, fun = max, na.rm = TRUE)

# bind all layers and extract all topographic variables at point location
topo_2m <- c(dem_vrt, slope_2m)
names(topo_2m) <- c("dem_2m","slope_2m")
vals_2m <- terra::extract(topo_2m, model_data)
model_data <- cbind(model_data, vals_2m[,-1])

topo_100m <- c(slope_100, aspect_100, dem_mean_100, dem_max_100, slope_max_100)
names(topo_100m) <- c("slope_100m", "aspect_100m","dem_mean_100m","dem_max_100m", "slope_max_100m")
vals_100m <- terra::extract(topo_100m, model_data)
model_data <- cbind(model_data, vals_100m[,-1])

# calculate orographic uplift at 100 m
model_data$OroUplift_100m <- oroUplift(u = model_data$U_10M_median,
                                       v = model_data$V_10M_median,
                                       slope_rad = (model_data$slope_100m * pi / 180), # in rad
                                       aspect_rad = (model_data$aspect_100m * pi / 180)) # in rad between 0 and 2pi
model_data$windDir_from_100m <- windDir_from_deg(u = model_data$U_10M_median,
                                                 v = model_data$V_10M_median)

ggplot(model_data, aes(x = aspect_100m, y = windDir_from_100m, color = OroUplift_100m)) +
  geom_point(size=1) +
  scale_color_gradient2(midpoint = 0) +
  coord_equal() +
  labs(x = "Slope aspect (degrees)",y = "Wind direction FROM (degrees)",fill = "w_oro (m/s)") +
  theme_minimal()

delta <- ((model_data$windDir_from_100m - model_data$aspect_100m + 180) %% 360) - 180
ggplot(model_data, aes(x = delta, y = OroUplift_100m)) +
  geom_point(alpha = 0.05) +
  geom_smooth() +
  theme_minimal() +
  labs(x = "Wind direction - aspect angular difference (deg)")

# calculate vertical distance of each nest presence/absence to the median elevation of the territory
table(model_data$territory_ID)
model_data <- model_data %>%
  group_by(territory_ID) %>%
  mutate(
    territory_median_elev = median(dem_2m, na.rm = TRUE),
    territory_max_elev    = max(dem_2m, na.rm = TRUE),
    # Differences from nest altitude relative to median and max altitude of territory
    verticalDist_nest_fromMedTerrElev = dem_2m - territory_median_elev,
    verticalDist_nest_fromMaxTerrElev    = dem_2m - territory_max_elev
  ) %>%
  ungroup()
class(model_data) <- setdiff(class(model_data), c("tbl_df", "tbl"))


# save model_data with the added variables
saveRDS(model_data, file="outputData/PCA_wind_pred_2708cio_allCRS_lulc_topo.rds")


