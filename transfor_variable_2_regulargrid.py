# Salish Sea Model Data Interpolation using Universal Kriging
# Introduction
# This python code  transform the output from the Salish Sea Model hydro model into
# regular grid (ROMS type) data using Universal Kriging. The Universal Kriging
# algorithm is employed to interpolate irregularly spaced data from the original
# Salish Sea Model grid to a regular grid, making it suitable for current tools
##
# **Author**: Javier Porobic
# **Date**: 01-07-2023
# **Last update**: 25-08-2023
##
# Background :
# The Salish Sea Model provides high-resolution data for the Salish Sea region,
# including temperature, salinity, and sigma layer values. The model data is often
# distributed on an irregular grid, which may not be directly suitable for specific
# modeling tasks or visualizations. To address this issue, we perform universal
# kriging to interpolate the data onto a regular grid.

# Dependencies
# - numpy
# - xarray
# - netCDF4
# - datetime
# - os
# - pyinterp

# Algorithm
# The main steps of the interpolation process are:
# 1. Read the original Salish Sea Model data from a NetCDF file.
# 2. Create a regular grid to interpolate the data onto.
# 3. Perform universal kriging to interpolate the temperature, salinity, and sigma layer values from the irregular grid to the regular grid.
# 4. Save the interpolated data into a new NetCDF file.


import numpy as np
import xarray as xr
import netCDF4
import datetime as dt
import os
import sys
import pyinterp
# Create an RTree instance for spatial indexing using pyinterp
mesh = pyinterp.RTree()

# Functions:
# Universal kriging interpolation


def kriging_universal(original_values, original_lon, original_lat, new_lon, new_lat):
    # Pack the original data into the RTree for spatial indexing
    mesh.packing(np.vstack((original_lon, original_lat)).T, original_values)
    # Perform universal kriging to interpolate new_lon, new_lat points
    kriging, neighbors = mesh.universal_kriging(np.vstack((new_lat.ravel(), new_lon.ravel())).T, within=True, k=11,
                                                covariance='matern_12', alpha=1_000_000, num_threads=0)
    return kriging.reshape(new_lon.shape)


# Define the file name of the original Salish Sea Model NetCDF data
filename = sys.argv[0]
filename_output = sys.argv[1]

print('Reading'+filename)
# Open the NetCDF file in read-write mode
# is necessary to change the name of the variable that is the same as the dimension
# this step is done just ones
# with netCDF4.Dataset(filename, mode='r+') as ds:
# Rename the variable that conflicts with a dimension name
#    ds.renameVariable('siglay', 'siglay_matrix')

# Open the NetCDF file and read the original Salish Sea Model data
ssm_solution = xr.open_dataset(filename, decode_cf=True, decode_times=False)
print('NetCDF file read!')
# MAIN GRID
# ~~~~~~~~~~~~~~
# Define the regular grid with a step of 0.01 for the new latitude and longitude
# (around 1km)
min_lat = ssm_solution.lat.min().values
max_lat = ssm_solution.lat.max().values
min_lon = ssm_solution.lon.min().values
max_lon = ssm_solution.lon.max().values
STEP = 0.01
reg_lat = np.arange(min_lon - STEP, max_lon + STEP, STEP)
reg_lon = np.arange(min_lat - STEP, max_lat + STEP, STEP)
# Meshgrid for the regular grid
mx, my = np.meshgrid(reg_lat, reg_lon, indexing='ij')
# Global Values
original_siglay = ssm_solution.siglay.values
original_siglev = ssm_solution.siglev.values
original_time = ssm_solution.time_vector.values
# Define the dimensions of the data
siglay_size = len(original_siglay)
siglev_size = len(original_siglev)
time_size = len(original_time)

# Variables
# Extract the original latitude, longitude, sigma layer, and time values
original_lat = ssm_solution.lat.values
original_lon = ssm_solution.lon.values


# Create empty arrays to store the interpolated temperature, salinity, and sigma layer values
new_regular_temp = np.full((len(original_time), len(
    original_siglay), len(reg_lat), len(reg_lon)), np.nan)
new_regular_salt = np.full((len(original_time), len(
    original_siglay), len(reg_lat), len(reg_lon)), np.nan)
new_regular_sigmalay = np.full(
    (len(original_siglay), len(reg_lat), len(reg_lon)), np.nan)

# Loop over each depth layer and interpolate the data onto the regular grid
for d in range(0, siglay_size):
    # Extract sigma layer values
    # siglev_matrix(siglev, node) ;
    org_sigma = ssm_solution.siglay_matrix[d].values
    new_regular_sigmalay[d][:] = kriging_universal(
        org_sigma, original_lon, original_lat, my, mx)

    for t in range(0, time_size):  # Loop over time steps
        org_temp = ssm_solution.temp[t][d].values  # Extract temperature values
        # Extract salinity values
        org_salt = ssm_solution.salinity[t][d].values
        new_regular_temp[t][d][:] = kriging_universal(
            org_temp, original_lon, original_lat, my, mx)
        new_regular_salt[t][d][:] = kriging_universal(
            org_salt, original_lon, original_lat, my, mx)

print('Interpolation variables done!')
# Velocity Field
# ~~~~~~~~~~~~~~
# Extract the original latitude, longitude, from each nele
original_latc = ssm_solution.latc.values
original_lonc = ssm_solution.lonc.values

# Create empty arrays to store the interpolated velocity fields
new_regular_v = np.full((len(original_time), len(
    original_siglay), len(reg_lat), len(reg_lon)), np.nan)
new_regular_u = np.full((len(original_time), len(
    original_siglay), len(reg_lat), len(reg_lon)), np.nan)

# Loop over each depth layer and interpolate the data onto the regular grid
for d in range(0, siglay_size):
    for t in range(0, time_size):  # Loop over time steps
        org_u = ssm_solution.u[t][d].values  # Extract temperature values
        org_v = ssm_solution.v[t][d].values  # Extract salinity values
        new_regular_u[t][d][:] = kriging_universal(
            org_u, original_lonc, original_latc, my, mx)
        new_regular_v[t][d][:] = kriging_universal(
            org_v, original_lonc, original_latc, my, mx)

print('Interpolation velocity field done!')

# Create a new NetCDF file with the interpolated data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save the interpolated data to a new NetCDF file
# data_path = r'/home/por07g/Documents/Projects/Salish_sea_puget_sound'
# file_name_output = 'regular_grid_variables_' + filename[-12:]
# nc = netCDF4.Dataset(os.path.join(data_path, file_name_output), 'w')
nc = netCDF4.Dataset(filename_output, 'w')

# Define NetCDF global attributes
nc.title = 'Regular grid for the Salish Sea model'
nc.Conventions = 'CF-1.6'
nc.history = '{0} creation of regular grid NetCDF file by Javier Porobic.'.format(
    dt.datetime.now().strftime("%Y-%m-%d"))

# Create NetCDF variables and set attributes
lat_dim = nc.createDimension('latitude', len(reg_lat))
lon_dim = nc.createDimension('longitude', len(reg_lon))
siglev_dim = nc.createDimension('sigma_layer', siglay_size)
time_dim = nc.createDimension('time', time_size)

lat_var = nc.createVariable('latitude', np.single, ('latitude'))
lat_var.units = 'degrees_north'
lat_var.standard_name = 'latitude'
lat_var.axis = 'Y'
lat_var[:] = reg_lat.astype('float')

lon_var = nc.createVariable('longitude', np.single, ('longitude'))
lon_var.units = 'degrees_east'
lon_var.standard_name = 'longitude'
lon_var.axis = 'X'
lon_var[:] = reg_lon.astype('float')

time_var = nc.createVariable('time_vector', np.intc, ('time'))
time_var.units = 'days since 1858-11-17 00:00:00'
time_var.standard_name = 'time'
time_var.format = 'modified julian day (MJD)'
time_var.time_zone = 'UTC'
time_var[:] = original_time.astype('int')

siglay_var = nc.createVariable('siglay', np.single, ('sigma_layer'))
siglay_var.units = 'sigma_layers'
siglay_var.standard_name = 'ocean_sigma/general_coordinate'
siglay_var[:] = original_siglay.astype('float')

siglev_var = nc.createVariable(
    'siglev', np.single, ('sigma_layer', 'latitude', 'longitude'))
siglev_var.units = 'sigma_level'
siglev_var.standard_name = 'ocean_sigma/general_coordinate'
siglev_var[:] = new_regular_sigmalay.astype('float')

temp_var = nc.createVariable(
    'temp', np.single, ('time', 'sigma_layer', 'latitude', 'longitude'))
temp_var.units = 'degrees_C'
temp_var.standard_name = 'sea_water_temperature'
temp_var[:] = new_regular_temp.astype('float')

salt_var = nc.createVariable(
    'salinity', np.single, ('time', 'sigma_layer', 'latitude', 'longitude'))
salt_var.units = '1e-3'
salt_var.standard_name = 'sea_water_salinity'
salt_var[:] = new_regular_salt.astype('float')

u_var = nc.createVariable(
    'u', np.single, ('time', 'sigma_layer', 'latitude', 'longitude'))
u_var.units = 'm s-1'
u_var.standard_name = 'eastward_sea_water_velocity'
u_var[:] = new_regular_u.astype('float')

v_var = nc.createVariable(
    'v', np.single, ('time', 'sigma_layer', 'latitude', 'longitude'))
v_var.units = 'm s-1'
v_var.standard_name = 'northward_sea_water_velocity'
v_var[:] = new_regular_v.astype('float')

nc.close()
print('NetCDF file created!')
