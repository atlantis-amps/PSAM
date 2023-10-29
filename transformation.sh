
## Time the variables in the file To simplified the output
ncks  -v temp,salinity,time,siglay,lon,lat SSM_v43_2018_0001.nc Variables.nc
## re-arranging the variables
ncpdq --rdr time,siglay,node Variables.nc -v  temp,salinity -o Variables2.nc
## removing the nodes
ncks -C -x -v node Variables2.nc -o clean_Variables.nc
## Interpolation between grids
cdo -L -remapnn,regular_grid.nc -selname,temp clean_Variables.nc nn_output.nc

#cdo remapbil,r360x180X10 SSM_v43_2018_0001.nc example.nc
cdo intlevel3d,r360x180X10 clean_Variables.nc example.nc
#Valocity
#ncks -v u,v SSM_v43_2018_0001.nc Velocity.nc



## help can be found here
##https://code.mpimet.mpg.de/boards/1/topics/14403?r=14741
