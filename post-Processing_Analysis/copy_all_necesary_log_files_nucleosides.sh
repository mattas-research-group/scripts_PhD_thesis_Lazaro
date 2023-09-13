#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid" "peptide")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a RU_names=("adenine" "guanine" "cytosine" "thymine" "uracil" "TARC_cbond" "TARC_nbond" "BA_cbond" "BA_nbond" "CA" "melamine")
declare -a RU_ac_names=("adenine_ac" "guanine_ac" "cytosine_ac" "thymine_ac" "uracil_ac" "TARC_cbond_ac" "TARC_nbond_ac" "BA_cbond_ac" "BA_nbond_ac" "CA_ac" "melamine_ac")


##########################################################################################################################
############################################## WORK ON EACH FOLDER #######################################################
##########################################################################################################################

# Open corresponding folder
cd nucleosides/final_gaussian_opt

for environment in "${environment[@]}"; do

cd $environment

for TC in "${TC_names[@]}"; do

cd $TC

if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

for ring_conf_fura in "${Ring_Conf_furanose[@]}"; do

cd $ring_conf_fura

for RU in "${RU_names[@]}"; do

cd $RU

# Copy TC final opt log file
cp $MD_RUN/nucleosides/initial_gaussian_opt/$environment/$TC/$ring_conf_fura/$RU/"$TC"_"$ring_conf_fura".log .
cp $MD_RUN/nucleosides/initial_gaussian_opt/$environment/$TC/$ring_conf_fura/$RU/"$RU".log .







cd ../
done

cd ../
done

elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then

for ring_conf_pyra in "${Ring_Conf_pyranose[@]}"; do

cd $ring_conf_pyra

for RU in "${RU_names[@]}"; do

cd $RU

# Copy TC final opt log file
cp $MD_RUN/nucleosides/initial_gaussian_opt/$environment/$TC/$ring_conf_pyra/$RU/"$TC"_"$ring_conf_pyra".log .
cp $MD_RUN/nucleosides/initial_gaussian_opt/$environment/$TC/$ring_conf_pyra/$RU/"$RU".log .









cd ../
done

cd ../
done

elif [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then

for RU in "${RU_names[@]}"; do

cd $RU

# Copy TC final opt log file
cp $MD_RUN/nucleosides/initial_gaussian_opt/$environment/$TC/$RU/"$TC".log .
cp $MD_RUN/nucleosides/initial_gaussian_opt/$environment/$TC/$RU/"$RU".log .









cd ../
done

elif [ "$TC" == "peptide" ]; then

for RU_ac in "${RU_ac_names[@]}"; do

cd $RU_ac

# Copy TC final opt log file
cp $MD_RUN/nucleosides/initial_gaussian_opt/$environment/$TC/$RU_ac/"$TC".log .
cp $MD_RUN/nucleosides/initial_gaussian_opt/$environment/$TC/$RU_ac/"$RU_ac".log .








cd ../
done

fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT
