#---------------------------------------------------------
# Import Modules
#---------------------------------------------------------
import sys
sys.path.append('/home/tcarrard/eagle.gravity.waves/utilities/')

import pandas as pd
import numpy as np

import xarray as xr
import xesmf as xe

from time import time

from rotate_grid import rotate_points

from multiprocessing import Pool

t0 = time()
#%%
#--------------------------------------------------------
# Settings
#--------------------------------------------------------

#define single levels interpolation data
interp_1l = ['U_10M','V_10M','ASHFL_S','HSURF','SLO_ASP_x', 'SLO_ASP_y', 'SLO_ANG_x', 'SLO_ANG_y']
interp_2l = ['U','V','W','TKE']

#--------------------------------------------------------
# Define time steps
#--------------------------------------------------------  

# Define start and end dates
start_date = pd.Timestamp("2023-05-01")
end_date = pd.Timestamp("2023-09-30")

# Generate base dates every 3 days
base_dates = pd.date_range(start=start_date, end=end_date, freq="3D")

# Define the times you want for each day
hours = [7, 13, 19]

# Generate full list of datetimes
datetime_list = [date + pd.Timedelta(hours=h) for date in base_dates for h in hours]

# Format as strings
formatted_times = [dt.strftime("%Y%m%d%H") for dt in datetime_list]

#--------------------------------------------------------
# Read data
#--------------------------------------------------------

eagle_path= '/net/thermo/atmosdyn/tcarrard/data/eagles/data/'
constpath   = '/net/litho/atmosdyn2/jansingl/KENDA-1/ANA23/'
constfile   = 'laf2023010100-const.nc'

filename = 'ForGiulia_backgandpresence_forAnnotation_Jul2025_noTime'

#read constant file
const_ds = xr.open_dataset(constpath+constfile)
z = const_ds['HEIGHT'].values
pole_lon = const_ds['grid_mapping_1'].grid_north_pole_longitude
pole_lat = const_ds['grid_mapping_1'].grid_north_pole_latitude

#define which constant variable need to be interpolated
interp_vars = const_ds[['HSURF','SLO_ASP','SLO_ANG']]
interp_vars[['SLO_ASP','SLO_ANG']] = interp_vars[['SLO_ASP','SLO_ANG']] * 180/np.pi

eagle_points = pd.read_csv(eagle_path+filename+'.csv')

#--------------------------------------------------------
# Subset COSMO
#--------------------------------------------------------

#subset COSMO in domain
min_lat = np.floor(eagle_points['y'].min())
max_lat =  np.ceil(eagle_points['y'].max())
min_lon =  np.floor(eagle_points['x'].min())
max_lon =  np.ceil(eagle_points['x'].max())

#rotate coordinates
lim_rlon,lim_rlat = rotate_points(pole_lon,pole_lat,[min_lon,max_lon],[min_lat,max_lat], 'n2r','NH')

#subset cosmo

cosmo_sub = const_ds.sel({'x_1':const_ds['x_1'].where((const_ds['x_1']>=lim_rlon[0]) & (const_ds['x_1'] <= lim_rlon[1]), drop=True)})
cosmo_sub1= cosmo_sub.sel({'y_1':cosmo_sub['y_1'].where((cosmo_sub['y_1']>=lim_rlat[0]) & (cosmo_sub['y_1'] <= lim_rlat[1]), drop=True)})

cosmo_sqz = cosmo_sub1.squeeze()

const_ds.close()

#--------------------------------------------------------
# Prepare grid for interpolation
#--------------------------------------------------------

eagle_points['rlon'], eagle_points['rlat'] = rotate_points(pole_lon,pole_lat,eagle_points['x'],eagle_points['y'],'n2r','NH')

cosmo_sqz = cosmo_sqz.rename({'x_1':'lon', 'y_1':'lat'})

chunk_size = 6000
num_chunks = int(np.ceil(len(eagle_points) / chunk_size))

regridders = []
eagle_subs = []

for i in range(num_chunks):
    batch = eagle_points.iloc[i * chunk_size : (i + 1) * chunk_size]

    target_grid = xr.Dataset({
        "lon": (("point1",), batch["rlon"].values),
        "lat": (("point2",), batch["rlat"].values)
    })
    
    print('Regridding ',chunk_size, ' points')

    regridder = xe.Regridder(
        cosmo_sqz,
        target_grid,
        method="bilinear",
        periodic=False,
        reuse_weights=False
    )

    regridders.append(regridder)
    eagle_subs.append(batch)
    
#--------------------------------------------------------
# Interpolating function
#--------------------------------------------------------  
step =formatted_times[0]

def interp_to_grid(step):

    #read cosmo file 
    
    cosmo_ds = xr.open_dataset(constpath+'laf'+step+'.nc')
    
    cosmo_ds.update(interp_vars)
    
    cosmo= cosmo_ds.sel({'x_1':cosmo_ds['x_1'].where((cosmo_ds['x_1']>=lim_rlon[0]) & (cosmo_ds['x_1'] <= lim_rlon[1]), drop=True)})
    cosmo1= cosmo.sel({'y_1':cosmo['y_1'].where((cosmo['y_1']>=lim_rlat[0]) & (cosmo['y_1'] <= lim_rlat[1]), drop=True)})
    
    cosmo_destaggered = cosmo1.rolling({'z_2': 2}, center=False).mean().isel({'z_2': slice(1, None)})
    
    #--------------------------------------------------------
    # Convert to single levels
    #--------------------------------------------------------
    
    #find which model level is closest to 200m
    hbs = z[:,:,:]-cosmo_ds['HSURF'].values
    hbs_lev = (hbs[:-1,:,:] + hbs[1:,:,:])/2
    hbs_lev = np.mean(hbs_lev,axis=(1,2)) #the model level height ranges between 150 and 250 meter above ground -> this is ok to take the average for our use.
    
    cosmo_ds.close()

    diff_height = abs(np.array([x-200 for x in hbs_lev]))
    
    lev_200 = np.argmin(diff_height)

    #subset cosmo levels between first level and lev 200
    cosmo2 = cosmo_destaggered.sel({'z_1':cosmo_destaggered['z_1'].where((cosmo_destaggered['z_1']>=lev_200), drop=True),
                                        'z_2':cosmo_destaggered['z_2'].where((cosmo_destaggered['z_2']>=lev_200), drop=True)})
    
    #average cosmo over model levels
    cosmo2[interp_2l] = cosmo2[interp_2l].mean(dim=['z_1','z_2'])
    
    #squeeze dataset 
    cosmo_squeezed = cosmo2.squeeze()
    
    cosmo_squeezed = cosmo_squeezed.rename({'x_1':'lon', 'y_1':'lat'})

    #--------------------------------------------------------
    # Switch slope angle to cartesian coordinates
    #--------------------------------------------------------
    
    #convert to cartesian coordinates
    cosmo_squeezed['SLO_ASP_x'] = np.cos(np.deg2rad(90-cosmo_squeezed['SLO_ASP']))
    cosmo_squeezed['SLO_ASP_y'] = np.sin(np.deg2rad(90-cosmo_squeezed['SLO_ASP']))
    cosmo_squeezed['SLO_ANG_x'] = np.cos(np.deg2rad(cosmo_squeezed['SLO_ANG']))
    cosmo_squeezed['SLO_ANG_y'] = np.sin(np.deg2rad(cosmo_squeezed['SLO_ANG']))
    
    for i,regridder in enumerate(regridders):
        
        print(round(i/len(regridders)*100,ndigits=1), 'completed')
        
        interpolated = regridder(cosmo_squeezed[interp_1l+interp_2l])        
        
        interpolated['SLO_ANG'] =  np.arctan2(interpolated['SLO_ANG_y'],interpolated['SLO_ANG_x']) #trigonometric in radians
        interpolated['SLO_ASP'] =  90-np.rad2deg(np.arctan2(interpolated['SLO_ASP_y'],interpolated['SLO_ASP_x'])) #cardinal in degrees
        
        #interpolated.drop_vars(['SLO_ANG_x','SLO_ANG_y','SLO_ASP_x','SLO_ASP_y'])
    
        #store as column in dataframe
        
        for var in interpolated.var():
            
            eagle_subs[i].loc[:,var] = interpolated[var].values.diagonal(axis1=-2, axis2=-1)
        
        #compute orographic lifting proxy
        eagle_subs[i]['UVspeed_10M'] = np.sqrt(eagle_subs[i]['U_10M']**2 + eagle_subs[i]['V_10M']**2)
        eagle_subs[i]['UVdir_10m'] = 270 - (180/np.pi)*np.arctan2(eagle_subs[i]['V_10M'],eagle_subs[i]['U_10M'])
        
        #calculate orographic uplift
        eagle_subs[i]['updraft_coeff'] = np.sin(eagle_subs[i]['SLO_ANG']) * np.cos(np.deg2rad(eagle_subs[i]['UVdir_10m'] - eagle_subs[i]['SLO_ASP']))
        eagle_subs[i]['w_oro'] = eagle_subs[i]['UVspeed_10M'] * eagle_subs[i]['updraft_coeff']
        eagle_subs[i]['time'] = step
        
    eagle_points = pd.concat(eagle_subs)

    return(eagle_points)

print('Starting interpolation for ', len(formatted_times), ' time steps')

#applly function to selected events
if __name__ == "__main__":
    with Pool(20) as pool:
      result = pool.map(interp_to_grid, formatted_times)
    print("Program finished!")

annotated_df = pd.concat(result)

print('Number of points annotated: ', annotated_df.index.size)

annotated_df.to_csv(eagle_path+filename[:-7]+'_annotated.csv')
    
tend = time()

print(tend-t0)

import pandas as pd
df = pd.read_csv('/home/tcarrard/eagle.gravity.waves/Scripts/ForGiulia_backgandpresence_forAnnotation_Jul2025_annotated.csv')

(df['UVdir_10m'] % 360).max()