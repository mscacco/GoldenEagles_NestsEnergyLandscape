#############################################################################
#  Script to produce a raster of distance from the edge of the territories  #
#############################################################################

#Loads libraries 
library(terra)
library(sf)
library(dplyr)
library(mapview)
library(Orcs)

#set my working directory 
setwd("F:/data")

sf_polygon <- readRDS("/data/tagged_territories.rds")
target_crs <- "EPSG:3035"
terr_sub3035 <-st_transform(sf_polygon, target_crs)
plot(terr_sub3035)

line_string <- st_cast(terr_sub3035, "MULTILINESTRING")
unique(line_string)
plot(line_string)

#load the base raster
dem_m <- rast("/data/raster_layer/final_rastercuttedongrigioni/dem_500m_resolution_3035.tif")
values(dem_m) <- NA


# Calculate distance from each cell in the raster to the nearest edge of the polygon (linestring)
distance_raster <- distance(dem_m, vect(line_string))
plot(distance_raster)
plot(line_string, add=T)


#save
writeRaster(distance_raster, filename = "/data/raster_layer/final_rastercuttedongrigioni/distance_to_edge_500mfinal.tif", overwrite=T)

check<- rast("/data/raster_layer/final_rastercuttedongrigioni/distance_to_edge_500mfinal.tif")
plot(check) #ok
