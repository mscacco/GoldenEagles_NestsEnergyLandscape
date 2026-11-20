#############################################################################
#  Script to produce a raster of distance from the edge of the territories  #
#############################################################################

#Loads libraries 
library(terra)
library(sf)
library(dplyr)
library(mapview)
library(Orcs)

#set working directory 
setwd("C:/Users/Giulia/Google Drive Streaming/My Drive/Paper nest selection Goldeneagle/GoldenEagles_NestsEnergyLandscape/Giulia_sharedFolder")

sf_nest_b <- read_sf("inputData/JuliaHatzl_nestdata/nestData_GR_onlypoints.gpkg") |>
  st_transform(3035)
sf_terr <- read_sf("inputData/JuliaHatzl_nestdata/territoriesGrisoni_Julia.gpkg") |>
  st_transform(3035)

# ---- Extract territory ID for each nest  --------------------------------------
hit   <- st_intersects(sf_nest, sf_terr)
first <- vapply(hit, function(ix) 
  if (length(ix)) ix[1] else NA_integer_, integer(1))
sf_nest$territory_id <- sf_terr$projectsite_id[first]   # index -> territory ID

terr_used <- sf_terr[sf_terr$projectsite_id %in% sf_nest$territory_id, ]

plot(terr_sub3035)

line_string <- st_cast(terr_used, "MULTILINESTRING")
unique(line_string)
plot(line_string)

#load the base raster
# Create a raster template 500 m EPSG:3035 to calculate distance ---------------
# Extent = territories + 5% on each side, snapped to multiples of res_m so the
# grid lines up with standard EPSG:3035 grids (e.g. the EEA reference grid).
res_m <- 500
e  <- ext(terr_used)
dx <- (xmax(e) - xmin(e)) * 0.05
dy <- (ymax(e) - ymin(e)) * 0.05
e_big <- ext(floor((xmin(e) - dx) / res_m) * res_m,
             ceiling((xmax(e) + dx) / res_m) * res_m,
             floor((ymin(e) - dy) / res_m) * res_m,
             ceiling((ymax(e) + dy) / res_m) * res_m)

r_template <- rast(e_big, resolution = res_m, crs = "EPSG:3035")
values(r_template) <- NA


# Calculate distance from each cell in the raster to the nearest edge of the polygon (linestring)
distance_raster <- distance(r_template, vect(line_string))
plot(distance_raster)
plot(line_string, add=T)

#save
writeRaster(distance_raster, filename = "outputData/distance_to_edge_500m.tif")

