netcdf regular_gird_variables_2018_0001 {
dimensions:
	latitude = 716 ;
	longitude = 785 ;
	sigma_layer = 10 ;
	time = 24 ;
variables:
	double latitude(latitude) ;
		latitude:units = "degrees_north" ;
		latitude:standard_name = "latitude" ;
		latitude:axis = "Y" ;
	double longitude(longitude) ;
		longitude:units = "degrees_east" ;
		longitude:standard_name = "longitude" ;
		longitude:axis = "X" ;
	double time_vector(time) ;
		time_vector:units = "days since 1858-11-17 00:00:00" ;
		time_vector:standard_name = "time" ;
		time_vector:format = "modified julian day (MJD)" ;
		time_vector:time_zone = "UTC" ;
	double siglay(sigma_layer) ;
		siglay:units = "sigma_layers" ;
		siglay:standard_name = "ocean_sigma/general_coordinate" ;
	double siglev(sigma_layer, latitude, longitude) ;
		siglev:units = "sigma_level" ;
		siglev:standard_name = "ocean_sigma/general_coordinate" ;
	double temp(time, sigma_layer, latitude, longitude) ;
		temp:units = "degrees_C" ;
		temp:standard_name = "sea_water_temperature" ;
	double salinity(time, sigma_layer, latitude, longitude) ;
		salinity:units = "1e-3" ;
		salinity:standard_name = "sea_water_salinity" ;

// global attributes:
		:title = "Regular grid for the Salish Sea model" ;
		:Conventions = "CF-1.6" ;
		:history = "2023-07-11 creation of regular grid netcdf file by Javier Porobic." ;
}
