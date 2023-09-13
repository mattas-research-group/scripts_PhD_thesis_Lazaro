#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a RU_names=("adenine" "guanine" "cytosine" "thymine" "uracil" "TARC_cbond" "TARC_nbond" "BA_cbond" "BA_nbond" "CA" "melamine")
declare -a IL_names=("phosphate" "arsenate")



function copy_files {

# Assign variables names to entries of function
local environment=$1
local TC=$2
local ring_conf=$3
local IL=$4
local RU=$5

if [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then 

# Copy files opt log for the RUs and the TCs
cp $MD_RUN/TCs_RUs_files/$environment/$TC/$RU/$RU.log .
cp $MD_RUN/TCs_RUs_files/$environment/$TC/$RU/$TC.log .

# Copy opt log files for the nucleoside and IL
cp $MD_RUN/gaussian_scan/$environment/$TC/$IL/$RU/"$TC"_"$RU".log .
# Copy opt log files for the IL
cp $MD_RUN/gaussian_scan/$environment/$TC/$IL/$RU/"$IL".log .

else

# Copy opt log files for the RUs and the TCs
cp $MD_RUN/TCs_RUs_files/$environment/$TC/$ring_conf/$RU/$RU.log .
cp $MD_RUN/TCs_RUs_files/$environment/$TC/$ring_conf/$RU/"$TC"_"$ring_conf".log .

# Copy opt log files for the nucleoside 
cp $MD_RUN/gaussian_scan/$environment/$TC/$ring_conf/$IL/$RU/"$TC"_"$ring_conf"_"$RU".log .
# Copy opt log files for the IL
cp $MD_RUN/gaussian_scan/$environment/$TC/$ring_conf/$IL/$RU/"$IL".log .

fi


}




##########################################################################################################################
############################################## WORK ON EACH FOLDER #######################################################
##########################################################################################################################

# Open corresponding folder
cd final_gaussian_opt

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

for RU in "${RU_names[@]}"; do

cd $RU

# Print message of on where folder is been worked
echo "working on $environment/$TC/$ring_conf_fura/$IL/$RU"

# Run script to obtain calculation status
copy_files $environment $TC $ring_conf_fura $IL $RU




cd ../
done

cd ../
done

cd ../
done

elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then

for ring_conf_pyra in "${Ring_Conf_pyranose[@]}"; do

cd $ring_conf_pyra

for IL in "${IL_names[@]}"; do

cd $IL

for RU in "${RU_names[@]}"; do

cd $RU

# Print message of on where folder is been worked
echo "working on $environment/$TC/$ring_conf_pyra/$IL/$RU"

# Run script to obtain calculation status
copy_files $environment $TC $ring_conf_pyra $IL $RU




cd ../
done

cd ../
done

cd ../
done

elif [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then

for IL in "${IL_names[@]}"; do

cd $IL

for RU in "${RU_names[@]}"; do

cd $RU

# Print message of on where folder is been worked
echo "working on $environment/$TC/$IL/$RU"

ring_conf="none"

# Run script to obtain calculation status
copy_files $environment $TC $ring_conf $IL $RU




cd ../
done

cd ../
done

fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT
