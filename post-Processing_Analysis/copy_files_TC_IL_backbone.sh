#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a IL_names=("phosphate" "arsenate")


############################################## WORK ON EACH FOLDER #######################################################

# Open corresponding folder
cd gaussian_scan

for environment in "${environment[@]}"; do

mkdir $environment

cd $environment

for TC in "${TC_names[@]}"; do

cd $TC

if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

for ring_conf_fura in "${Ring_Conf_furanose[@]}"; do

cd $ring_conf_fura

for IL in "${IL_names[@]}"; do

cd $IL

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, ring conf $ring_conf_fura and IL $IL"

# Copy files
cp $IL.log $MD_RUN/final_gaussian_opt/$environment/$TC/$ring_conf_fura/$IL/$IL.log
cp "$TC"_"$ring_conf_fura".log $MD_RUN/final_gaussian_opt/$environment/$TC/$ring_conf_fura/$IL/"$TC"_"$ring_conf_fura".log






cd ../
done

cd ../
done

elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then

for ring_conf_pyra in "${Ring_Conf_pyranose[@]}"; do

cd $ring_conf_pyra

for IL in "${IL_names[@]}"; do

cd $IL


# Print message of on where folder is been worked
echo "working on $environment, TC $TC, ring conf $ring_conf_pyra and IL $IL"

# Copy files
cp $IL.log $MD_RUN/final_gaussian_opt/$environment/$TC/$ring_conf_pyra/$IL/$IL.log
cp "$TC"_"$ring_conf_pyra".log $MD_RUN/final_gaussian_opt/$environment/$TC/$ring_conf_pyra/$IL/"$TC"_"$ring_conf_pyra".log
 




cd ../
done

cd ../
done

elif [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then

for IL in "${IL_names[@]}"; do

cd $IL

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, none ring conf and IL $IL"

# Copy files
cp $IL.log $MD_RUN/final_gaussian_opt/$environment/$TC/$IL/$IL.log
cp "$TC".log $MD_RUN/final_gaussian_opt/$environment/$TC/$IL/"$TC".log







cd ../
done


fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT
