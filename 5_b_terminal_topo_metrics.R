library(terra)
library(terra)

# ---- settings ----

dem_fld <- "C:/Users/Giulia/Desktop/topo_work/DEM_Swisstopo_Grisons_2m/"
out_dir <- "C:/Users/Giulia/Desktop/topo_work/outputs/"
tmp_dir <- "C:/Users/Giulia/Desktop/topo_work/tmp_terra/"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

terraOptions(tempdir = tmp_dir, memfrac = 0.3, progress = 1)
wopt <- list(datatype = "FLT4S",
             gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))

f <- function(name) file.path(out_dir, paste0(name, ".tif"))

# run a step only if its output doesn't exist yet, then free memory/temp files
step <- function(name, fn) {
  if (file.exists(f(name))) { message("skip  ", name); return(rast(f(name))) }
  message(format(Sys.time()), "  start ", name)
  r <- fn(f(name))
  message(format(Sys.time()), "  done  ", name)
  gc()
  r
}

# ---- input ----
tiles <- list.files(dem_fld, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
message("Number of tiles found: ", length(tiles))
if (length(tiles) == 0) stop("No .tif tiles found in ", dem_fld)
dem_vrt <- vrt(tiles, filename = file.path(out_dir, "dem_2m.vrt"), overwrite = TRUE)

# ---- 2 m slope (the big one; everything fine-scale uses it) ----
slope_2m <- step("slope_2m", \(fn) terrain(dem_vrt, v = "slope", unit = "degrees",
                                           filename = fn, wopt = wopt))

# ---- 100 m ----
dem_mean_100  <- step("dem_mean_100",  \(fn) aggregate(dem_vrt, 50, fun = "mean", na.rm = TRUE, filename = fn, wopt = wopt))
slope_100     <- step("slope_100",     \(fn) terrain(dem_mean_100, v = "slope",  unit = "degrees", filename = fn, wopt = wopt))
aspect_100    <- step("aspect_100",    \(fn) terrain(dem_mean_100, v = "aspect", unit = "degrees", filename = fn, wopt = wopt))
slope_max_100 <- step("slope_max_100", \(fn) aggregate(slope_2m, 50, fun = "max", na.rm = TRUE, filename = fn, wopt = wopt))
dem_max_100   <- step("dem_max_100",   \(fn) aggregate(dem_vrt,  50, fun = "max", na.rm = TRUE, filename = fn, wopt = wopt))

# ---- 10 m ----
dem_mean_10  <- step("dem_mean_10",  \(fn) aggregate(dem_vrt, 5, fun = "mean", na.rm = TRUE, filename = fn, wopt = wopt))
slope_10     <- step("slope_10",     \(fn) terrain(dem_mean_10, v = "slope",  unit = "degrees", filename = fn, wopt = wopt))
aspect_10    <- step("aspect_10",    \(fn) terrain(dem_mean_10, v = "aspect", unit = "degrees", filename = fn, wopt = wopt))
slope_max_10 <- step("slope_max_10", \(fn) aggregate(slope_2m, 5, fun = "max", na.rm = TRUE, filename = fn, wopt = wopt))
dem_max_10   <- step("dem_max_10",   \(fn) aggregate(dem_vrt,  5, fun = "max", na.rm = TRUE, filename = fn, wopt = wopt))

message("ALL DONE ", format(Sys.time()))