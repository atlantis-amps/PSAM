#!/bin/bash
## Activate your environment in my case is GBR_env
conda activate GBR_env
# Set the input and output folder names. This will need to be changed for your system
root_folder="/home/por07g/Documents/Projects/Salish_sea_puget_sound"
input_folder="$root_folder/input_Data"
copy_folder="$root_folder/copy_Data"
output_folder="$root_folder/interpolated_data"
## creating the output and copy folder
mkdir -p "$output_folder"
mkdir -p "$copy_folder"
# Loop over the input files in the input folder
for input_file in "$input_folder"/*.nc; do
    # Get the base name of the input file (without the path or extension)
    input_basename=$(basename "$input_file" .nc)
    echo "Reading... $input_basename"
    # Set the output file name based on the input file name
    output_file="$output_folder/${input_basename}_interpolated.nc"
    # Make a copy of the input file with a new name
    cp "$input_file" "$copy_folder/${input_basename}_copy.nc"
    # # create the new input file name    
    new_input_file="$copy_folder/${input_basename}_copy.nc"
    echo $new_input_file
    # Run the Python script with the input and output file names as arguments
    # you should use the updated python script that we workend on the other day
    python transfor_variable_2_regulargrid.py "$new_input_file" "$output_file"
done

## these are the commands that I used to create the average files or the sampling every 12 hours

# ncra -y avg -d time,0,11 -d time,12,23 "$output_folder/SSM_v43_2018_0001_original.nc" "$output_folder/SSM_v43_2018_0001_avg.nc"
## ncks -d time,0,,2 input_file.nc extracted_file.nc
