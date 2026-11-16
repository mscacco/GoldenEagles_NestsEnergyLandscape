#_____________________________________________________________#
#Annotate orographic uplift resolution at xx m resolution   ###
#_____________________________________________________________#


#Script to recalculate orographic uplift at higher resolution

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
