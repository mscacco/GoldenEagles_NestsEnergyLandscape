#___________________________________________________
# Annotate topographic and land cover variables ####
#___________________________________________________

# This script derives and aggregates topographic layers from 2 m pixels to 100 m and extract summary information to highlight the presence of cliffs and surface structures
# It annotates both topographic, land cover, and distance from territory data to the dataset
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


setwd("/home/martina/ownCloud/Martina/ProgettiVari/GoldenEagles_COSMO/Nests_EnergyLandscape_Giulia/Giulia_sharedFolder")
setwd("C:/Users/Giulia/Google Drive Streaming/My Drive/Paper nest selection Goldeneagle/GoldenEagles_NestsEnergyLandscape/Giulia_sharedFolder")

# Load dataset (sf object) with all CRS
model_data <- readRDS("outputData/df_all_coordinates_forAnnotation.rds")


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
writeRaster(lulc_stack_100m, "inputData/environmental_data/LandCover_10m_PlanComp/lulc_derived100m_pcts.tif", overwrite=TRUE)

# annotate data with info from all layers
model_data_prj <- st_transform(model_data, crs(lucl_10m))
val_lulc_10m <- terra::extract(lucl_10m, model_data_prj, method="nearest")
val_lulc_100m <- terra::extract(lulc_stack_100m, model_data_prj, method="nearest")

model_data_lulc <- cbind(model_data, lulc_10m=val_lulc_10m[,-1], val_lulc_100m[,-1])
table(model_data_lulc$lulc_10m)
table(model_data_lulc$lulc_100m) # they are very similar, we retain the one at 100m


# add legend and save
length(unique(model_data_lulc$ID_Tom))==nrow(model_data_lulc) #we have to check this 
model_data_lulc <- merge(model_data_lulc, leg, by.x="lulc_100m", by.y="value")
model_data_lulc <- dplyr::rename(model_data_lulc, lcClass_100m = class)
table(model_data_lulc$lulc_100m)
table(model_data_lulc$lcClass_100m)

#save
saveRDS(model_data_lulc, "outputData/df_allCRS_lulc.rds") #before was PCA_wind_pred_2708cio_allCRS_lulc.rds


#____________________________________
# Derive and annotate topography ####

# import the data annotated with the land cover
model_data <- readRDS("outputData/df_allCRS_lulc.rds") #before was ( PCA_wind_pred_2708cio_allCRS_lulc.rds)

## load all tiles from HD to R as a virtual raster (does NOT load all into memory)
##output_fld <- "/media/martina/Intenso/DEM_Swisstopo_Grisons_2m/"
#output_fld <- "inputData/environmental_data/DEM_Swisstopo_Grisons_2m/"
#tiles <- list.files(output_fld, pattern = "\\.tif", full.names = TRUE)
#dem_vrt <- terra::vrt(tiles)

## Calculate 4 complementary topographic metrics (100m): 
## Broad-scale slope (1) and aspect (2), calculated from a 100 m DEM obtained by aggregating the original 2 m DEM using the mean
#dem_mean_100 <- terra::aggregate(dem_vrt, fact = 50, fun = mean, na.rm = TRUE)
#slope_100 <- terra::terrain(dem_mean_100, v = "slope", unit = "degrees")
#aspect_100 <- terra::terrain(dem_mean_100, v = "aspect", unit = "degrees")
#
## Fine-scale terrain ruggdeness, calculated as the maximum dem (3) and maximum slope (4) within 100 m windows based on dem at 2 m and slope computed at 2 m resolution.
#slope_2m <- terra::terrain(dem_vrt, v = "slope", unit = "degrees")
#slope_max_100 <- terra::aggregate(slope_2m, fact = 50, fun = max, na.rm = TRUE)
#dem_max_100 <- terra::aggregate(dem_vrt, fact = 50, fun = max, na.rm = TRUE)
#
## Calculate 4 complementary topographic metrics (10m): 
## Broad-scale slope (1) and aspect (2), calculated from a 10 m DEM obtained by aggregating the original 2 m DEM using the mean
#dem_mean_10 <- terra::aggregate(dem_vrt, fact = 5, fun = mean, na.rm = TRUE)
#slope_10 <- terra::terrain(dem_mean_10, v = "slope", unit = "degrees")
#aspect_10 <- terra::terrain(dem_mean_10, v = "aspect", unit = "degrees")
#
## Fine-scale terrain ruggdeness, calculated as the maximum dem (3) and maximum slope (4) within 10 m windows based on dem at 2 m and slope computed at 2 m resolution.
#slope_max_10 <- terra::aggregate(slope_2m, fact = 5, fun = max, na.rm = TRUE)
#dem_max_10 <- terra::aggregate(dem_vrt, fact = 5, fun = max, na.rm = TRUE)
#
## bind all layers and extract all topographic variables at point location
#
##2m
#topo_2m <- c(dem_vrt, slope_2m)
#names(topo_2m) <- c("dem_2m","slope_2m")
#vals_2m <- terra::extract(topo_2m, model_data)
#model_data <- cbind(model_data, vals_2m[,-1])
#
##10m
#topo_10m <- c(slope_10,aspect_10, dem_mean_10, dem_max_10, slope_max_10)
#names(topo_10m) <- c("slope_10m", "aspect_10m", "dem_mean_10m", "dem_max_10m", "slope_max_10m")
#vals_10m <- terra::extract(topo_10m, model_data)
#
##100m
#topo_100m <- c(slope_100, aspect_100, dem_mean_100, dem_max_100, slope_max_100)
#names(topo_100m) <- c("slope_100m", "aspect_100m","dem_mean_100m","dem_max_100m", "slope_max_100m")
#vals_100m <- terra::extract(topo_100m, model_data)
#model_data <- cbind(model_data, vals_100m[,-1])

#____________________________________

#if you run the code on terminal using code topo_metrics.R run the following 

library(terra)
library(sf)
library(dplyr)

setwd("C:/Users/Giulia/Google Drive Streaming/My Drive/Paper nest selection Goldeneagle/GoldenEagles_NestsEnergyLandscape/Giulia_sharedFolder")

# import the data annotated with the land cover
model_data <- readRDS("outputData/df_allCRS_lulc.rds")

# topographic rasters computed overnight (topo_metrics.R)
topo_dir <- "C:/Users/Giulia/Desktop/topo_work/outputs/"
ld <- function(name) rast(file.path(topo_dir, paste0(name, ".tif")))

# 2 m DEM (rebuilt from the local tiles, takes a few seconds)
tiles   <- list.files("C:/Users/Giulia/Desktop/topo_work/DEM_Swisstopo_Grisons_2m/",
                      pattern = "\\.tif$", full.names = TRUE)
dem_vrt  <- vrt(tiles)
slope_2m <- ld("slope_2m")

dem_mean_100 <- ld("dem_mean_100"); slope_100 <- ld("slope_100"); aspect_100 <- ld("aspect_100")
slope_max_100 <- ld("slope_max_100"); dem_max_100 <- ld("dem_max_100")

dem_mean_10 <- ld("dem_mean_10"); slope_10 <- ld("slope_10"); aspect_10 <- ld("aspect_10")
slope_max_10 <- ld("slope_max_10"); dem_max_10 <- ld("dem_max_10")

#points in the same CRS as the DEM before extracting
model_data_prj <- st_transform(model_data, crs(dem_vrt))

#_________________________
# Extract topographic variables at point locations and add them to model_data

#2m
topo_2m <- c(dem_vrt, slope_2m)
names(topo_2m) <- c("dem_2m", "slope_2m")
model_data <- cbind(model_data, terra::extract(topo_2m, model_data_prj, ID = FALSE))

#10m
topo_10m <- c(slope_10, aspect_10, dem_mean_10, dem_max_10, slope_max_10)
names(topo_10m) <- c("slope_10m", "aspect_10m", "dem_mean_10m", "dem_max_10m", "slope_max_10m")
model_data <- cbind(model_data, terra::extract(topo_10m, model_data_prj, ID = FALSE))

#100m
topo_100m <- c(slope_100, aspect_100, dem_mean_100, dem_max_100, slope_max_100)
names(topo_100m) <- c("slope_100m", "aspect_100m", "dem_mean_100m", "dem_max_100m", "slope_max_100m")
model_data <- cbind(model_data, terra::extract(topo_100m, model_data_prj, ID = FALSE))

#Should we calculate aspect in this way?
# aspect as northness / eastness (aspect is circular)
model_data$northness_10m  <- cos(model_data$aspect_10m  * pi / 180)
model_data$eastness_10m   <- sin(model_data$aspect_10m  * pi / 180)
model_data$northness_100m <- cos(model_data$aspect_100m * pi / 180)
model_data$eastness_100m  <- sin(model_data$aspect_100m * pi / 180)



# quick checks
c("territory_id", "dem_2m") %in% names(model_data)
summary(model_data[, c("dem_2m", "slope_10m", "slope_100m")])

#where do we have NAs
has_na <- !complete.cases(st_drop_geometry(model_data)[, c("dem_2m", "slope_10m", "slope_100m")])
table(case = model_data$case, has_na = has_na)

footprint  <- as.polygons(ifel(is.na(dem_vrt), NA, 1))
na_pts <- vect(model_data[model_data$has_na, ])

plot(dem_vrt)
plot(footprint, add = TRUE, border = "black", lwd = 1.5, col = NA)
plot(na_pts, add = TRUE, pch = 21, bg = "red", col = "black", cex = 0.2, lwd = 0.3)

#_________________________
# calculate vertical distance of each nest presence/absence to the median elevation of the territory
table(model_data$territory_id, useNA = "ifany")
model_data <- model_data %>%
  group_by(territory_id) %>%
  mutate(
    territory_median_elev = median(dem_2m, na.rm = TRUE),
    territory_max_elev    = max(dem_2m, na.rm = TRUE),
    # Differences from nest altitude relative to median and max altitude of territory
    verticalDist_nest_fromMedTerrElev = dem_2m - territory_median_elev,
    verticalDist_nest_fromMaxTerrElev = dem_2m - territory_max_elev
  ) %>%
  ungroup()
class(model_data) <- setdiff(class(model_data), c("tbl_df", "tbl"))


#_________________________
#Distance from territory boundaries 

dem_dist <- rast("outputData/distance_to_edge_500m.tif")
names(dem_dist) <- "dist_terr_500m"
pts <- st_transform(model_data, crs(dem_dist))
model_data <- cbind(model_data, terra::extract(dem_dist, pts, ID = FALSE))


# save model_data with the added variables
saveRDS(model_data, file="outputData/model_data_allCRS_lulc_topo_dist.rds")



