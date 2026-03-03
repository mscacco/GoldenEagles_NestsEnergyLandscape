#______________________________________________#
#Computing mean variables annotated with COSMO #
#______________________________________________#


#Script to compute mean value of the variables annotated using COSMO 
#we will have mean of the following variables:
#- W vertical wind speed 
#- UVdir horizontal wind direction
#- UVspeed horizontal wind speed 
#- oroUplift orographic uplift 
#- TKE turbulence kinetic energy 
#- surface sensible heat flux



#load dataset annotated by Tom


#remove HSURF, slope, aspect since are calcualted at resolution of 1.1km


#compute mean 


#add covariates to previous dataframe 