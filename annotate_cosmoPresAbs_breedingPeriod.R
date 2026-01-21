
#_______________________________
# Prepare COSMO annotation ####
#_______________________________

# This script prepares the dataset for annotation with atmospheric variables from COSMO
# While our dataset is static, the atmospheric variables are dynamic in time
# We were interested in the variability, predictability and average atmospheric conditions at the nest location (presences) relative to absence locations during the breeding season
# The script associates the same presence and absence location with a sequence of timesteps that represent the breeding season
# To reduce computation time it takes 3 hours per day every third day across the season (May to September)
# The output dataset is then annotated using the python script "interpolate_average_cosmo.py" at ETH, where COSMO data are locally stored


library(lubridate)
library(sf)
library(mapview)
library(data.table)

setwd("/home/martina/ownCloud/Martina/ProgettiVari/GoldenEagles_COSMO/Nests_EnergyLandscape_Giulia/Giulia_sharedFolder/")


#_______________________
# Define timesteps ####

# Define the breeding period
breedingPeriod <- as.Date(c("01-05-2023","30-09-2023"), format="%d-%m-%Y")

# Make a sequence of all days within this period
allDays <- seq(breedingPeriod[1], breedingPeriod[2], by=3) # Every third day, to reduce computation time
# Define hours of interest (sunrise to sunset, three time points per day to reduce computation time)
allTimes <- paste0(sprintf("%02d", c(7,13,19)), ":00") #make sure all have 2 digits (leading 0s)

# All combinations of dates and times as POSIXct object
allTimestamps <- as.POSIXct(apply(expand.grid(allDays, allTimes), 1, paste, collapse=" "),
                            format="%Y-%m-%d %H:%M", tz="CET")
allTimestamps <- with_tz(allTimestamps, "UTC")
length(allTimestamps)
range(allTimestamps)

#___________________________________________
# Associate annotation times to dataset ####

# Import presence absence
presAbs <- readRDS("inputData/final_df_23072025_quota_forAnnotation.rds")

# Project to lat long decimal wgs84
presAbs_wgs84 <- st_transform(presAbs, st_crs("EPSG:4326"))
mapview(presAbs_wgs84)

presAbs_wgs84$x <- st_coordinates(presAbs_wgs84$geometry)[,1]
presAbs_wgs84$y <- st_coordinates(presAbs_wgs84$geometry)[,2]
PAdf <- st_drop_geometry(presAbs_wgs84)
PAdf <- as.data.frame(PAdf[order(PAdf$territory_id, -PAdf$case),])
head(PAdf)
table(PAdf$case)

# Explore:
# N of presences per territory (between 1 and 14 eagles' nests per territory)
table(PAdf[PAdf$case==1,"territory_id"])
table(table(PAdf[PAdf$case==1,"territory_id"]))
# N of territories (121)
length(unique(PAdf$territory_id))
# N of points per territory (between 70 and > 700 points)
range(sapply(split(PAdf, PAdf$territory_id), nrow))

# Repeat points for each timestamp during the breeding period
nrow(PAdf) * length(allTimestamps) # expected number of rows of final dataset

PAdf_time <- as.data.frame(rbindlist(lapply(allTimestamps, function(time){
  return(cbind(PAdf, timestamp = time))
})))
nrow(PAdf_time) == nrow(PAdf) * length(allTimestamps) 
length(unique(PAdf_time$timestamp))
range(PAdf_time$timestamp)
unique(table(PAdf_time$timestamp))

write.csv(PAdf_time, file = "outputData/finaldataset_backgandpresence_timeSteps_forAnnotation_Jul2025.csv", row.names = F)
saveRDS(PAdf_time, file = "outputData/finaldataset_backgandpresence_timeSteps_forAnnotation_Jul2025.rds")

# the final annotated dataset will be named "ForGiulia_backgandpresence_timeSteps_forAnnotation_Jul2025_annotated.csv"
