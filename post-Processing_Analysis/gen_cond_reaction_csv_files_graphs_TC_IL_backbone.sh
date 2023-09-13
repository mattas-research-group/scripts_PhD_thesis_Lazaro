#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a IL_names=("phosphate" "arsenate")


function num_gaussian_opt {

# Assign variables names to entries of function
local environment=$1
local TC=$2
local ring_conf=$3
local IL=$4

if [ "$ring_conf" == "none" ]; then
log_file_prefix="$TC"_"$IL"
else
log_file_prefix="$TC"_"$ring_conf"_"$IL"
fi

# Run if then statement to define maximum number of opt+freq
find -type f -name "$log_file_prefix*" -and -name "*min.log" > list_gaussian_min.txt
find -type f -name "$log_file_prefix*" -and -name "*min2.log" > list_gaussian_opt.txt
sort -n list_gaussian_opt.txt > tmp.txt && mv tmp.txt list_gaussian_opt.txt
cat list_gaussian_min.txt list_gaussian_opt.txt > final_list_gaussian.txt

# Delete first two characters ./ from each line
cut -c 3- final_list_gaussian.txt > tmp.txt && mv tmp.txt list_gaussian.txt

# Select last file to work on it
tail -n 1 list_gaussian.txt > final_list_gaussian.txt
rm list_gaussian.txt

}


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

function calc_thermo_params_cond_reaction {

# Assign variables names to entries of function
local log_file=$1
local component_name=$2
local component_type=$3
local TC=$4
local destfolder=$5

########################################################## Find thermo properties #######################################################
# Print number of line of last instance of Maximum Force
awk '/Sum of electronic and zero-point Energies=/ { ln = FNR } END { print ln }' "$log_file" > line_num_zeropoint_energy.txt

############################################# OBTAIN ZERO POINT ENERGY ################################

# Obtain the line contents 
line_n_zeropoint_energy=`cat line_num_zeropoint_energy.txt`
sed -n ''$line_n_zeropoint_energy'p' "$log_file" > zeropoint_energy.txt

# Print 2nd column with actual energy value
awk '{ print $7 }' zeropoint_energy.txt > zeropoint_energy_value.txt

# Convert energy to kcal/mol
awk '{ printf "%.8e\n", $1*2625.5 }' zeropoint_energy_value.txt > zeropoint_energy_kj_mol.txt

# Make copy of this file with appropiate name for component
cp zeropoint_energy_kj_mol.txt "$component_name"_zeropoint_energy_kj_mol.txt

# Assign variable to zero point energy in kcal/mol txt file
zero_point_energy=`cat "$component_name"_zeropoint_energy_kj_mol.txt`

############################################# OBTAIN CORRECTED ENERGY ################################
line_corrected_energy=$(($line_n_zeropoint_energy+1))
sed -n ''$line_corrected_energy'p' "$log_file" > corrected_energy_react.txt

# Print 2nd column with actual energy value
awk '{ print $7 }' corrected_energy_react.txt > corrected_energy_value_react.txt

# Convert energy to kcal/mol
awk '{ printf "%.8e\n", $1*2625.5 }' corrected_energy_value_react.txt > corrected_energy_kj_mol.txt

# Make copy of this file with appropiate name for component
cp corrected_energy_kj_mol.txt "$component_name"_corrected_energy_kj_mol.txt

# Assign variable to zero point energy in kcal/mol txt file
corrected_energy=`cat "$component_name"_corrected_energy_kj_mol.txt`

############################################# OBTAIN ENTHALPY ################################
line_enthalpy=$(($line_n_zeropoint_energy+2))
sed -n ''$line_enthalpy'p' "$log_file" > entalphy_react.txt

# Print 2nd column with actual energy value
awk '{ print $7 }' entalphy_react.txt > enthalpy_value_react.txt

# Convert energy to kcal/mol
awk '{ printf "%.8e\n", $1*2625.5 }' enthalpy_value_react.txt > entalphy_kj_mol.txt

# Make copy of this file with appropiate name for component
cp entalphy_kj_mol.txt "$component_name"_entalphy_kj_mol.txt

# Assign variable to zero point energy in kcal/mol txt file
enthalpy=`cat "$component_name"_entalphy_kj_mol.txt`

############################################ OBTAIN THE FREE ENERGY ##################################
line_dG=$(($line_n_zeropoint_energy+3))
sed -n ''$line_dG'p' "$log_file" > dG_react.txt

# Print 2nd column with actual energy value
awk '{ print $8 }' dG_react.txt > dG_value_react.txt

# Convert energy to kcal/mol
awk '{ printf "%.8e\n", $1*2625.5 }' dG_value_react.txt > G_kj_mol.txt

# Make copy of this file with appropiate name for component
cp G_kj_mol.txt "$component_name"_G_kj_mol.txt

# Assign variable to zero point energy in kcal/mol txt file
free_energy=`cat "$component_name"_G_kj_mol.txt`

# OUTPUT ENERGIES AND OTHER INFO IN CORRESPONDING SUMMARY FILES IN ANALYSIS FOLDER
if [ "$component_type" == "backbone" ]; then
echo $zero_point_energy >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/backbone_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/backbone_corrected_energies.txt
echo $enthalpy >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/backbone_enthalpies.txt
echo $free_energy >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/backbone_free_energies.txt

elif [ "$component_type" == "IL" ]; then
echo $zero_point_energy >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/ILs_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/ILs_corrected_energies.txt
echo $enthalpy >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/ILs_enthalpies.txt
echo $free_energy >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/ILs_free_energies.txt

elif [ "$component_type" == "TC" ]; then
echo $zero_point_energy >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/TCs_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/TCs_corrected_energies.txt
echo $enthalpy >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/TCs_enthalpies.txt
echo $free_energy >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/TCs_free_energies.txt
fi


}


function output_other_info {

# Declare variables
local environment=$1
local TC=$2
local ring_conf=$3
local IL=$4
local destfolder=$5

# Output names of components in general files
echo $environment >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/environment_list.txt
echo $TC >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/TCs_names_list.txt
echo $ring_conf >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/ring_conf_names_list.txt
echo $IL >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/ILs_names_list.txt


}




function output_results_for_water {

local environment=$1
local destfolder=$2

# Water parameters
if [ "$environment" == "vacuum" ]; then
zpE_water=$(echo "scale=8; -76.437249*2625.5" | bc)
corrE_water=$(echo "scale=8; -76.434413*2625.5" | bc)
H_water=$(echo "scale=8; -76.433469*2625.5" | bc)
G_water=$(echo "scale=8; -76.454892*2625.5" | bc)
elif [ "$environment" == "water" ]; then
zpE_water=$(echo "scale=8; -76.445286*2625.5" | bc)
corrE_water=$(echo "scale=8; -76.442450*2625.5" | bc)
H_water=$(echo "scale=8; -76.441506*2625.5" | bc)
G_water=$(echo "scale=8; -76.462936*2625.5" | bc)
fi

echo $zpE_water >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/water_zero_point_energies.txt
echo $corrE_water >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/water_corrected_energies.txt
echo $H_water >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/water_enthalpies.txt
echo $G_water >> $MD_RUN/post_proc_summary_results_csv_graphs/condensation_reaction/$destfolder/water_free_energies.txt

}







############################################## WORK ON EACH FOLDER #######################################################

# Create folder for analysis and post-processing
mkdir post_proc_summary_results_csv_graphs
mkdir post_proc_summary_results_csv_graphs/condensation_reaction
mkdir post_proc_summary_results_csv_graphs/condensation_reaction/furanoses
mkdir post_proc_summary_results_csv_graphs/condensation_reaction/pyranoses
mkdir post_proc_summary_results_csv_graphs/condensation_reaction/non_sugar_TC

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

# Run function to determine how many opt+freq were runned and get log file TC_IL_backbone
num_gaussian_opt $environment $TC $ring_conf_fura $IL 

# Get variable for final log file of backbone
log_file=`cat final_list_gaussian.txt`
destfolder="furanoses"


# Run function to extract thermo parameters from log files for TC_IL_backbone, IL and TC 
calc_thermo_params_cond_reaction $log_file "$TC"_"$ring_conf_fura"_"$IL" "backbone" $TC $destfolder
calc_thermo_params_cond_reaction $IL.log $IL "IL" $TC $destfolder
calc_thermo_params_cond_reaction "$TC"_"$ring_conf_fura".log "$TC"_"$ring_conf_fura" "TC" $TC $destfolder

# Output name of components
output_other_info $environment $TC $ring_conf_fura $IL $destfolder

# Output energies from water
output_results_for_water $environment $destfolder



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
 
# Run function to determine how many opt+freq were runned and get log file for TC_IL_backbone
num_gaussian_opt $environment $TC $ring_conf_pyra $IL 

# Get variable for final log file of backbone
log_file=`cat final_list_gaussian.txt`
destfolder="pyranoses"

# Run function to extract thermo parameters from log files for TC_IL_backbone, IL and TC 
calc_thermo_params_cond_reaction $log_file "$TC"_"$ring_conf_pyra"_"$IL" "backbone" $TC $destfolder
calc_thermo_params_cond_reaction $IL.log $IL "IL" $TC $destfolder
calc_thermo_params_cond_reaction "$TC"_"$ring_conf_pyra".log "$TC"_"$ring_conf_pyra" "TC" $TC $destfolder

# Output name of components
output_other_info $environment $TC $ring_conf_pyra $IL $destfolder

# Output energies from water
output_results_for_water $environment $destfolder







cd ../
done

cd ../
done

elif [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then

for IL in "${IL_names[@]}"; do

cd $IL

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, none ring conf and IL $IL"

ring_conf="none"

# Run function to determine how many opt+freq were runned and get log file of TC_IL backbone
num_gaussian_opt $environment $TC $ring_conf $IL 

# Get variable for final log file of backbone
log_file=`cat final_list_gaussian.txt`
destfolder="non_sugar_TC"

# Run function to extract thermo parameters from log files for TC_IL_backbone, IL and TC 
calc_thermo_params_cond_reaction $log_file "$TC"_"$IL" "backbone" $TC $destfolder
calc_thermo_params_cond_reaction $IL.log $IL "IL" $TC $destfolder
calc_thermo_params_cond_reaction "$TC".log "$TC" "TC" $TC $destfolder

# Output name of components
ring_conf="none"
output_other_info $environment $TC $ring_conf $IL $destfolder

# Output energies from water
output_results_for_water $environment $destfolder



cd ../
done


fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT
