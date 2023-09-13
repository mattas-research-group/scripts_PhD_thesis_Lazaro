#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")


##########################################################################################################################
############################################## WORK ON EACH FOLDER #######################################################
##########################################################################################################################

# Open corresponding folder
mkdir TC
cd TC
mkdir final_gaussian_opt
cd final_gaussian_opt

for environment in "${environment[@]}"; do

mkdir $environment
cd $environment

for TC in "${TC_names[@]}"; do

mkdir $TC
cd $TC

if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

for ring_conf_fura in "${Ring_Conf_furanose[@]}"; do

mkdir $ring_conf_fura
cd $ring_conf_fura

# Copy TC final opt log file
cp $MD_RUN/nucleosides/initial_gaussian_opt/$environment/$TC/$ring_conf_fura/adenine/"$TC"_"$ring_conf_fura".log .


cd ../ # LEAVE CURRENT RING CONF

done

else

for ring_conf_pyra in "${Ring_Conf_pyranose[@]}"; do

mkdir $ring_conf_pyra
cd $ring_conf_pyra

# Copy TC final opt log file
cp $MD_RUN/nucleosides/initial_gaussian_opt/$environment/$TC/$ring_conf_pyra/adenine/"$TC"_"$ring_conf_pyra".log .



cd ../ # LEAVE CURRENT RING CONF

done

fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT
