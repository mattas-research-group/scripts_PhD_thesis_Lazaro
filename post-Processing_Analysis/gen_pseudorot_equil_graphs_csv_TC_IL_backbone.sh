#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a IL_names=("phosphate" "arsenate")


##################################################
##                                              ##
## FUNCTIONS TO OBTAIN DIFFERENT THERMODYNAMIC  ##
## PARAMETERS FOR CONDENSATION REACTION.        ##
##                                              ##
## Description:                                 ##
## This function extract E, H and G from log    ##
## files of nucleoside, TC and RU and put them  ##
## in corresponding nucleoside file.            ##
##                                              ##
## It also create summary txt files in analysis ##
## folder.                                      ##
##                                              ##
##################################################

function extract_thermo_params_pseu_rot {

# Assign variables names to entries of function
local environment=$1
local TC=$2
local IL=$3
local ring_conf=$4
local destfolder=$5
local component_name=$6

############################################# OBTAIN ZERO POINT ENERGY ################################

# Assign variable to zero point energy in kcal/mol txt file
zero_point_energy=`cat "$component_name"_zeropoint_energy_kj_mol.txt`

############################################# OBTAIN CORRECTED ENERGY ################################

# Assign variable to zero point energy in kcal/mol txt file
corrected_energy=`cat "$component_name"_corrected_energy_kj_mol.txt`

############################################# OBTAIN ENTHALPY ################################

# Assign variable to zero point energy in kcal/mol txt file
enthalpy=`cat "$component_name"_entalphy_kj_mol.txt`

############################################ OBTAIN THE FREE ENERGY ##################################

# Assign variable to zero point energy in kcal/mol txt file
free_energy=`cat "$component_name"_G_kj_mol.txt`

# OUTPUT ENERGIES AND OTHER INFO IN CORRESPONDING SUMMARY FILES IN ANALYSIS FOLDER
echo $environment >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$destfolder/"$ring_conf"_environment_list.txt
echo $TC >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$destfolder/"$ring_conf"_TCs_names_list.txt
echo $IL >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$destfolder/"$ring_conf"_ILs_names.txt
echo $ring_conf >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$destfolder/"$ring_conf"_ring_conf_names_list.txt

echo $zero_point_energy >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$destfolder/"$ring_conf"_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$destfolder/"$ring_conf"_corrected_energies.txt
echo $enthalpy >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$destfolder/"$ring_conf"_enthalpies.txt
echo $free_energy >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$destfolder/"$ring_conf"_free_energies.txt


}





############################################## WORK ON EACH FOLDER #######################################################

# Create folder for analysis and post-processing
mkdir post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium
mkdir post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/furanoses
mkdir post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/pyranoses

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

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, ring conf $ring_conf_fura and IL $IL" 

# Get variable for final log file of backbone
log_file=`cat final_list_gaussian.txt`
destfolder="furanoses"

# Run function to extract thermo parameters from log files for nucleoside, RU and TC 
extract_thermo_params_pseu_rot $environment $TC $IL $ring_conf_fura $destfolder "$TC"_"$ring_conf_fura"_"$IL"



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

# Get variable for final log file of backbone
log_file=`cat final_list_gaussian.txt`
destfolder="pyranoses"

# Run function to extract thermo parameters from log files for nucleoside, RU and TC 
extract_thermo_params_pseu_rot $environment $TC $IL $ring_conf_pyra $destfolder "$TC"_"$ring_conf_pyra"_"$IL"







cd ../
done

cd ../
done

fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT
