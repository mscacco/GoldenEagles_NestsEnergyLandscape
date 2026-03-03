
#__________________________________________
# Download high-res DEM and land cover ####
#__________________________________________

# This script downloads a 2 m DEM for the study area (Grisons, Switzerland)
# https://www.swisstopo.admin.ch/en/height-model-swissalti3d
# And a 10 m land cover from the planetary computer
# https://planetarycomputer.microsoft.com/dataset/io-lulc-annual-v02
# Aggregates the pixels and extracts summary information to highlight the presence of cliffs and surface structures


library(rstac)
library(terra)
library(sf)
library(future)
library(future.apply)
library(curl)

setwd("/home/martina/ownCloud/Martina/ProgettiVari/GoldenEagles_COSMO/Nests_EnergyLandscape_Giulia/Giulia_sharedFolder/")

# Load dataset (sf object)
model_data <- readRDS("outputData/PCA_wind_pred_2708cio.rds")

head(model_data)
names(model_data)

#_________________________________
# CRSs and BBOXs for download ####

# store laea coordinates
crsLaea <- crs(model_data)
model_data$x_laea <- st_coordinates(model_data)[,1]
model_data$y_laea <- st_coordinates(model_data)[,2]

# Reproject to wgs84 and store coordinates and extent
model_data <- st_transform(model_data, 4326)
crswgs84 <- crs(model_data)
model_data$x_wgs84 <- st_coordinates(model_data)[,1]
model_data$y_wgs84 <- st_coordinates(model_data)[,2]

(bbox_wgs84 <- ext(model_data)*1.1)


# Reproject to swiss coordinate system LV95 (EPSG:2056)
model_data <- st_transform(model_data, 2056)
crsLv95 <- crs(model_data)
model_data$x_lv95 <- st_coordinates(model_data)[,1]
model_data$y_lv95 <- st_coordinates(model_data)[,2]

(bbox_lv95 <- ext(model_data)*1.1)

# resave model_data with new added coordinates and crs
saveRDS(model_data, "outputData/PCA_wind_pred_2708cio_allCRS.rds")


#__________________________________
# Download Swisstopo DEM tiles ####

# We used this bbox to manually download the list of urls of the swiss topo DEM at 2 m.
# These can be freely downloaded from:
# https://www.swisstopo.admin.ch/en/height-model-swissalti3d

bbox_lv95

# The urls of the >10000 1km tiles to download are stored in the inputData folder, but the tiles will be downloaded on an external hard drive
output_fld <- "/media/martina/Intenso/DEM_Swisstopo_Grisons_2m/"
urls <- read.csv("inputData/environmental_data/downloadLinks_tilesGrisons_swisstopoSwissalti3d_modelDataExt.csv")[,1]

# if only partially downloaded run this to continue download
table(basename(urls) %in% list.files(output_fld, pattern="tif"))
urls <- urls[!basename(urls) %in% list.files(output_fld, pattern="tif")]

downloadTile <- function(u, output_dir) {
  dest <- file.path(output_fld, basename(u))
  if (file.exists(dest)) return(dest)
  
  tryCatch({
    curl_download(u, destfile = dest, quiet = TRUE)
    dest
  }, error = function(e) {
    message("Failed: ", u)
    NA
  })
}

# Run in Terminal
plan(multicore, workers = future::availableCores() - 2)
files <- future_lapply(urls, downloadTile, output_fld)
# to check for download errors
# files <- unlist(files)
# files <- files[!is.na(files)]


#__________________________
# Download Land Cover ####

# We use the bbox in lat long to download the land cover from the planetary computer via rstac
bbox_wgs84 <- as.vector(bbox_wgs84)[c("xmin", "ymin", "xmax", "ymax")]

output_dir <- "inputData/environmental_data/LandCover_10m_PlanComp/"

stac_source <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")
collections_query <- stac_source |> collections()
available_collections <- get_request(collections_query)
collections_title <- sapply(available_collections$collections, function(x) x$title)
collections_title[grep("Land Use Land Cover", collections_title)]
collections_id <- sapply(available_collections$collections, function(x) x$id)
collections_id[collections_title == "10m Annual Land Use Land Cover (9-class) V2"]

# "10m Annual Land Use Land Cover (9-class) V2"
# https://planetarycomputer.microsoft.com/api/stac/v1/collections/io-lulc-annual-v02


stac_query <- stac_search(
  q = stac_source, # Stac query object
  collections = "io-lulc-annual-v02",
  datetime = "2023-01-01/2023-12-31",
  bbox = bbox_wgs84
)
(executed_stac_query <- get_request(stac_query))
filtered_stac_query <- items_select(executed_stac_query, 1) # select first feature (only year 2023)
items_properties(filtered_stac_query)

# Sign the stac query - this is necessary to actually download the data - 
signed_stac_query <- items_sign(
  filtered_stac_query,
  sign_planetary_computer()
)
signed_urls <- assets_url(signed_stac_query)

# It is only one tile, the first file is a raster and the second is the metadata
download.file(
  url = signed_urls[1],
  destfile = paste0(output_dir, "lulc_v02_2023.tif"),   
  mode = "wb"
)
download.file(
  url = signed_urls[3],
  destfile = paste0(output_dir, "metadata.json"),   
  mode = "wb"
)

r <- rast(paste0(output_dir, "lulc_v02_2023.tif"))
plot(r)

# legend available directly at https://planetarycomputer.microsoft.com/dataset/io-lulc-annual-v02
# and saved as txt file:
(leg <- read.table(paste0(output_dir, "lulc_v02_legend.txt"), sep=",", header = T))
