## in order to let CDO work with your input, the dimension of the data
## variable mi need to be re-ordered. this is possible with NCO:


# Next is the deletion of the data variable ncells:

ncks -C -x -v ncells output.nc -o cleaned_output.nc

# then you can interpolate with something like:

cdo -L -remapnn,global_1 -selname,mi cleaned_output.nc nn_output.nc

