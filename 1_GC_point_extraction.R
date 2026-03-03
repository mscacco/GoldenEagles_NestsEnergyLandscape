#__________________________________________
# Extract absence points ####
#__________________________________________


#This script load presence data and territory boundaries to extract absence points 
#
#
#Assign weight to points according to dimension of territories 



#load libraries 
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(units)
library(car)
library(caret)
library(patchwork)
library(factoextra)
library(mapview)
library(mgcv)
library(tidyr)
library(amt)

#set working directory 
setwd("F:/data") #setwd(.....percorso fino alla cartella condivisa)

#dalla cartella condivisa in avanti 

#load data and extract info about territories 
#brood_data <- read.csv("/script/golden_eagle_MSc_old/golden_eagle_MSc/Paper/New folder/brood_data2023_grisons.csv", sep = ";")
sf_nest <- read_sf("/script/golden_eagle_MSc_old/golden_eagle_MSc/Paper/New folder/nestData_GR_onlypoints.gpkg")  %>% st_transform("EPSG:3035")
sf_terr<- read_sf("territoriesGrisoni_Julia.gpkg") %>% st_transform("EPSG:3035")

##load my predictors and align them together 
#bands <- rast("F:/data/raster_layer/final_rastercuttedongrigioni/sentinel2_mean_500m_resolution_3035.tif")
#bands_1 <- bands[[1]]
#
#dtm <- rast("F:/data/raster_layer/dem_m.tif")
#rough <- rast("F:/data/raster_layer/roughness_30.tif")
#slope <- rast("F:/data/raster_layer/slope_30.tif")
#aspect <- rast("F:/data/raster_layer/aspect_30.tif") 
#
#distance <- rast("F:/data/raster_layer/final_rastercuttedongrigioni/distance_raster_aligned_500m_resolution_3035.tif")
#
#
##stack rasters together
#stack <- c(rough, slope, aspect,dtm)
##stack bands dem distance
#stack_2 <- c(distance, bands)
#crop_stack <- crop(stack, stack_2)
#names(crop_stack)[4] <- "dtm"


#######presence and bg extraction #################
#Rasterise territories and stack them with other layer since we need to extract absences from center of each cell that have not presence point  
# Convert the shapefile to a SpatVector (terra object)
terr_vect <- vect(sf_terr)

#extract terr info 
intersection_list <- st_intersects(sf_nest, sf_terr)

matched_ids <- sapply(intersection_list, function(ix) {
  if (length(ix) == 0) return("CHECK")  # No intersection
  return(sf_terr$projectsite_id[ix[1]])  # Use first match if multiple
})

# Add result to the points data
sf_nest$matched_polygon_id <- matched_ids


unique(terr_vect$projectsite_id)
ids_to_keep <- unique(sf_nest$matched_polygon_id) #terr ID to keep 

terr_vect_to_keep <- terr_vect[terr_vect$projectsite_id %in% ids_to_keep, ] #filter ID to keep from terr_vect
plot(terr_vect_to_keep)

#get the extension from sf_terr
extension <- ext(terr_vect_to_keep)

#create the raster template 
r_template <- bands[[1]]
values(r_template) <- NA #assign NA values 

#extract info about territories for the nest 
#check CRS they need to have the same 
st_crs(sf_nest) == st_crs(sf_terr) #ok


# Extract only centroids that are within the polygons

# r_template is the template raster with the same structure of the other env layers, poly is the sf with polygons
# Extract raster centroids
r_df <- xyFromCell(r_template, 1:ncell(r_template))
plot(r_template)
plot(r_df$geometry, add=T)

r_df_sf <- st_as_sf(as.data.frame(r_df), coords= c("x", "y"), crs= st_crs(r_template))

# Extract only centroids that are within the polygons
df_poly <- r_df_sf %>%
  st_join(st_as_sf(terr_vect_to_keep), join = st_within, left = FALSE)
# left = FALSE is an inner_join, so keeps only points that have a matching polygon

# Extract coordinates and add as new columns
df_poly <- df_poly %>%
  mutate(x = st_coordinates(geometry)[, 1],
         y = st_coordinates(geometry)[, 2]) %>% 
  select(-geometry) %>% 
  st_drop_geometry() 

# Convert to data.table (geometry is preserved unless you drop it)
df_poly <- as.data.table(df_poly) %>% 
  select(-name)

# rasterise next presence points on same r_template
presence_point <- sf_nest %>% 
  st_transform("EPSG:3035") %>% 
  select(geom, matched_polygon_id) %>% 
  rename(territory_id = matched_polygon_id) %>% 
  mutate(case = 1)

presence_r <- rasterize(vect(presence_point), r_template, touches = TRUE, background = NA)

# extract coordinates of pixels including presence points
non_na_cells <- which(!is.na(values(presence_r))) #cells with presence (non-NA)
presence_centroids <- xyFromCell(presence_r, non_na_cells) # get centroids of pixels containing nests
# Exclude those presence centroids from the list of background points
presence_centroids <- as.data.table(presence_centroids)
setnames(presence_centroids, c("x", "y"))  # match column names in df_poly
background_centroids <- anti_join(df_poly, presence_centroids, by = c("x", "y")
)

# Sanity check
nrow(df_poly)               # original number of background centroids
nrow(presence_centroids)    # number of presence points
nrow(background_centroids)   # after filtering

#add info about case to background_centroids
background_centroids$case <- 0

#add info about terr it to presence 
presence_point <- presence_point %>% 
  rename(geometry = geom)

df_final <- bind_rows(presence_point, background_point) #27909

#(df_final, "/script/golden_eagle_MSc_old/golden_eagle_MSc/Paper/New folder/final_df_2307.rds")
#table(df_final$territory_id[df_final$case == 1])
#

df_final_vect <- vect(df_final)

#extract the quota just in case it's possible to annotate at the nest level
df_quota <- terra::extract(dem, df_final_vect)

df_final_q <- bind_cols(df_final, df_quota)
saveRDS(df_final_q, "/script/golden_eagle_MSc_old/golden_eagle_MSc/Paper/New folder/final_df_2307_quota.rds")
#################################

#load df final and first part about the raster 
df <- readRDS("F:/script/golden_eagle_MSc_old/golden_eagle_MSc/Paper/New folder/final_df_2307_quota.rds")


