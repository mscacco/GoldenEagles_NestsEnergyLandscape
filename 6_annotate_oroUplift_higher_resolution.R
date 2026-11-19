#_____________________________________________________________#
#Annotate orographic uplift resolution at xx m resolution   ###
#_____________________________________________________________#


#Script to compute mean value of the variables annotated using COSMO and 
#recalculate orographic uplift at higher resolution 

#we will have mean of the following variables:
#- W vertical wind speed 
#- UVdir horizontal wind direction
#- UVspeed horizontal wind speed 
#- oroUplift orographic uplift 
#- TKE turbulence kinetic energy 
#- surface sensible heat flux


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

#load model data annotate with topo and land and merge with df from tom 

model_data <- readRDS("outputData/model_data_allCRS_lulc_topo_dist.rds")
cosmo_data <- read.csv()

#
#remove HSURF, slope, aspect since are calcualted at resolution of 1.1km
#compute mean, median ecc 


#join them together 
model_data <- model_data %>% 
  left_join(cosmo_aggregated_data, ecc )


# calculate orographic uplift at 100 m
model_data$OroUplift_100m <- oroUplift(u = model_data$U_10M_median,
                                       v = model_data$V_10M_median,
                                       slope_rad = (model_data$slope_100m * pi / 180), # in rad
                                       aspect_rad = (model_data$aspect_100m * pi / 180)) # in rad between 0 and 2pi

#calculate orographic uplift at 10 m
model_data$OroUplift_10m <- oroUplift(u = model_data$U_10M_median,
                                      v = model_data$V_10M_median,
                                      slope_rad = (model_data$slope_10m * pi / 180), # in rad
                                      aspect_rad = (model_data$aspect_10m * pi / 180)) # in rad between 0 and 2pi


#should we delete this? I guess is going to provite it 
#Wind direction (we can omit 100m since doesn't really have a resolution right? but just to keep in mind was recalculate from us an not tom)
model_data$windDir_from_100m <- windDir_from_deg(u = model_data$U_10M_median,
                                                 v = model_data$V_10M_median)
#___________________________


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
