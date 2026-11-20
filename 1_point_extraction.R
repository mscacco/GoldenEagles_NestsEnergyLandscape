############### Presence and background extraction (500 m grid) ###############
# Output: one row per presence (exact nest location, case = 1) and one row per
# background cell centre inside a used territory with no nest (case = 0),
# plus elevation from the DTM. Everything is in EPSG:3035 (metres).

library(sf)
library(terra)
library(data.table)

#set working directory 
setwd("C:/Users/Giulia/Google Drive Streaming/My Drive/Paper nest selection Goldeneagle/GoldenEagles_NestsEnergyLandscape/Giulia_sharedFolder")


res_m <- 500   # grid resolution in metres (old dataset used ~497.94 m)

# ---- 1. Load data -------------------------------------------------------------
sf_nest <- read_sf("inputData/JuliaHatzl_nestdata/nestData_GR_onlypoints.gpkg") |>
  st_transform(3035)
sf_terr <- read_sf("inputData/JuliaHatzl_nestdata/territoriesGrisoni_Julia.gpkg") |>
  st_transform(3035)

# ---- Extract territory ID for each nest  --------------------------------------
hit   <- st_intersects(sf_nest, sf_terr)
first <- vapply(hit, function(ix) 
  if (length(ix)) ix[1] else NA_integer_, integer(1))
sf_nest$territory_id <- sf_terr$projectsite_id[first]   # index -> territory ID

#check NAs
n_out <- sum(is.na(sf_nest$territory_id))
if (n_out > 0) {
  warning(n_out, " nest(s) fall outside every territory and are dropped. ",
          "Inspect them with: sf_nest[is.na(sf_nest$territory_id), ]")
  sf_nest <- sf_nest[!is.na(sf_nest$territory_id), ]
}

# Keep only territories that contain at least one nest
terr_v <- vect(sf_terr[sf_terr$projectsite_id %in% sf_nest$territory_id, ])

# Create a raster template 500 m EPSG:3035 to extract background coordinates ---------------
# Extent = territories + 5% on each side, snapped to multiples of res_m so the
# grid lines up with standard EPSG:3035 grids (e.g. the EEA reference grid).
e  <- ext(terr_v)
dx <- (xmax(e) - xmin(e)) * 0.05
dy <- (ymax(e) - ymin(e)) * 0.05
e_big <- ext(floor((xmin(e) - dx) / res_m) * res_m,
             ceiling((xmax(e) + dx) / res_m) * res_m,
             floor((ymin(e) - dy) / res_m) * res_m,
             ceiling((ymax(e) + dy) / res_m) * res_m)

r_template <- rast(e_big, resolution = res_m, crs = "EPSG:3035")
ncell(r_template)

# Extract coordinates fro each cells inside territory ------------------
# A cell in two overlapping territories appears once per territory.
cand <- as.data.table(cells(r_template, terr_v))
cand[, territory_id := terr_v$projectsite_id[ID]]   # index -> territory ID
cand[, ID := NULL]

# Combine Presences and background points ---------------------------------------------
nest_xy <- st_coordinates(sf_nest)

pres <- data.table(
  territory_id = sf_nest$territory_id,
  cell         = cellFromXY(r_template, nest_xy[, 1:2]),
  x            = nest_xy[, 1],
  y            = nest_xy[, 2],
  case         = 1L
)

bg <- cand[!cell %in% pres$cell]                     # drop cells holding a nest
bg[, c("x", "y") := as.data.table(xyFromCell(r_template, cell))]
bg[, case := 0L]

df_final <- rbind(pres, bg, use.names = TRUE)

# Sanity checks (old dataset: 602 presences, 27307 background)
cat("Candidate cells:  ", nrow(cand), "\n",
    "Unique cells:     ", uniqueN(cand$cell), "\n",
    "Presences:        ", nrow(pres), "\n",
    "Background cells: ", nrow(bg), "\n")
print(df_final[, .N, by = case])

# Visual check
plot(terr_v, border = "grey40")
points(bg$x, bg$y, pch = ".", col = "grey50")
points(pres$x, pres$y, pch = 16, cex = 0.5, col = "red")

# ---- 6. Elevation from DTM ---------------------------------------------------
dtm <- rast("inputData/environmental_data/dtm_elev_eumap_epsg3035_v0.3_cropAlps_laea.tif")
plot(dtm)
pts <- vect(df_final, geom = c("x", "y"), crs = "EPSG:3035")
df_final[, elev := extract(dtm, project(pts, crs(dtm)), ID = FALSE)[[1]]]

saveRDS(df_final, "outputData/final_df_500m_quota_forAnnotation.rds")
################################################################################