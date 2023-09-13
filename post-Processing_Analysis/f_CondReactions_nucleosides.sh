#!/bin/bash

#	f_CondReactions_nucleosides.sh

#	08-jun-2023	Edited options for subtitle and group files in TC folders in python function
# 	13-feb-2023	Created

# Description:
# This script calculates the values of all the thermo properties for all nucleosides and
# generates all graphs and summary csv files. 


# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid" "peptide")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a RU_names=("adenine" "guanine" "cytosine" "thymine" "uracil" "TARC_cbond" "TARC_nbond" "BA_cbond" "BA_nbond" "CA" "melamine")
declare -a RU_ac_names=("adenine_ac" "guanine_ac" "cytosine_ac" "thymine_ac" "uracil_ac" "TARC_cbond_ac" "TARC_nbond_ac" "BA_cbond_ac" "BA_nbond_ac" "CA_ac" "melamine_ac")


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

local log_file=$1
local component_name=$2
local component_type=$3


########################################################## Find thermo properties #######################################################
# Print number of line of last instance of Maximum Force
awk '/Sum of electronic and zero-point Energies=/ { ln = FNR } END { print ln }' "$log_file" > line_num_zeropoint_energy.txt

############################################# OBTAIN ZERO POINT ENERGY ################################

# Obtain the line contents 
line_n_zeropoint_energy=`cat line_num_zeropoint_energy.txt`
sed -n ''$line_n_zeropoint_energy'p' "$log_file" > zeropoint_energy.txt

# Print 2nd column with actual energy value
awk '{ print $7 }' zeropoint_energy.txt > zeropoint_energy_value.txt

# Convert energy to kJ/mol
awk '{ printf "%.8e\n", $1*627.503 }' zeropoint_energy_value.txt > zeropoint_energy_kJ_mol.txt

# Make copy of this file with appropiate name for component
cp zeropoint_energy_kJ_mol.txt "$component_name"_zeropoint_energy_kJ_mol.txt

# Assign variable to zero point energy in kJ/mol txt file
zero_point_energy=`cat "$component_name"_zeropoint_energy_kJ_mol.txt`

############################################# OBTAIN CORRECTED ENERGY ################################
line_corrected_energy=$(($line_n_zeropoint_energy+1))
sed -n ''$line_corrected_energy'p' "$log_file" > corrected_energy_react.txt

# Print 2nd column with actual energy value
awk '{ print $7 }' corrected_energy_react.txt > corrected_energy_value_react.txt

# Convert energy to kJ/mol
awk '{ printf "%.8e\n", $1*627.503 }' corrected_energy_value_react.txt > corrected_energy_kJ_mol.txt

# Make copy of this file with appropiate name for component
cp corrected_energy_kJ_mol.txt "$component_name"_corrected_energy_kJ_mol.txt

# Assign variable to zero point energy in kJ/mol txt file
corrected_energy=`cat "$component_name"_corrected_energy_kJ_mol.txt`

############################################# OBTAIN ENTHALPY ################################
line_enthalpy=$(($line_n_zeropoint_energy+2))
sed -n ''$line_enthalpy'p' "$log_file" > entalphy_react.txt

# Print 2nd column with actual energy value
awk '{ print $7 }' entalphy_react.txt > enthalpy_value_react.txt

# Convert energy to kJ/mol
awk '{ printf "%.8e\n", $1*627.503 }' enthalpy_value_react.txt > entalphy_kJ_mol.txt

# Make copy of this file with appropiate name for component
cp entalphy_kJ_mol.txt "$component_name"_entalphy_kJ_mol.txt

# Assign variable to zero point energy in kJ/mol txt file
enthalpy=`cat "$component_name"_entalphy_kJ_mol.txt`

############################################ OBTAIN THE FREE ENERGY ##################################
line_dG=$(($line_n_zeropoint_energy+3))
sed -n ''$line_dG'p' "$log_file" > dG_react.txt

# Print 2nd column with actual energy value
awk '{ print $8 }' dG_react.txt > dG_value_react.txt

# Convert energy to kJ/mol
awk '{ printf "%.8e\n", $1*627.503 }' dG_value_react.txt > G_kJ_mol.txt

# Make copy of this file with appropiate name for component
cp G_kJ_mol.txt "$component_name"_G_kJ_mol.txt

# Assign variable to zero point energy in kJ/mol txt file
free_energy=`cat "$component_name"_G_kJ_mol.txt`



# OUTPUT ENERGIES AND OTHER INFO IN CORRESPONDING SUMMARY FILES IN ANALYSIS FOLDER
if [ "$component_type" == "nucleoside" ]; then
echo $zero_point_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/nucleosides_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/nucleosides_corrected_energies.txt
echo $enthalpy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/nucleosides_enthalpies.txt
echo $free_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/nucleosides_free_energies.txt

elif [ "$component_type" == "TC" ]; then
echo $zero_point_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/TCs_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/TCs_corrected_energies.txt
echo $enthalpy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/TCs_enthalpies.txt
echo $free_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/TCs_free_energies.txt

else

echo $zero_point_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/RUs_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/RUs_corrected_energies.txt
echo $enthalpy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/RUs_enthalpies.txt
echo $free_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/RUs_free_energies.txt

fi

#################################################################################### END OF FUNCTION ###########################################################################################

}

function output_other_info {

local ring_conf=$1
local RU=$2

echo $RU >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/RUs_names_list.txt
echo $TC >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/TCs_names_list.txt
echo $ring_conf >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/ring_conf_names_list.txt
echo $environment >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/environment_list.txt

}


function output_results_for_water {

local environment=$1

# Water parameters
if [ "$environment" == "vacuum" ]; then
zpE_water=$(echo "scale=8; -76.437249*627.503" | bc)
corrE_water=$(echo "scale=8; -76.434413*627.503" | bc)
H_water=$(echo "scale=8; -76.433469*627.503" | bc)
G_water=$(echo "scale=8; -76.454892*627.503" | bc)
elif [ "$environment" == "water" ]; then
zpE_water=$(echo "scale=8; -76.445286*627.503" | bc)
corrE_water=$(echo "scale=8; -76.442450*627.503" | bc)
H_water=$(echo "scale=8; -76.441506*627.503" | bc)
G_water=$(echo "scale=8; -76.462936*627.503" | bc)
fi

echo $zpE_water >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/water_zero_point_energies.txt
echo $corrE_water >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/water_corrected_energies.txt
echo $H_water >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/water_enthalpies.txt
echo $G_water >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/water_free_energies.txt

}

##########################################################################################################################
############################################## WORK ON EACH FOLDER #######################################################
##########################################################################################################################

# Create folder for analysis and post-processing
mkdir nucleosides/post_proc_summary_results_csv_graphs
mkdir nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction

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

# Output message of status of script
echo "working on $environment/$TC/$ring_conf_fura/$RU"

# Define if there is one or two or three opt log files for nucleoside
# Define input file for function that calculates the thermo equilibrium
log_file_prefix="$TC"_"$ring_conf_fura"_"$RU"

# Define input file for function that calculates the thermo equilibrium
unique_opt_file="$log_file_prefix"_min.log
second_opt_file="$log_file_prefix"_2opt.log
third_opt_file="$log_file_prefix"_3opt.log

if [ -f "$third_opt_file" ]; then

log_file=$third_opt_file

echo "3opt" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/n_gaussian_opt.txt

else

if [ -f "$second_opt_file" ]; then

log_file=$second_opt_file

echo "2opt" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/n_gaussian_opt.txt

else

log_file=$unique_opt_file

echo "1opt" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/n_gaussian_opt.txt

fi

fi

# Run function to extract thermo parameters from log files for nucleoside, RU and TC 
component_TC="TC"
component_RU="RU"
component_nucleoside="nucleoside"
TC_name="$TC"_"$ring_conf_fura"
nucleoside_name="$TC"_"$ring_conf_fura"_"$RU"
calc_thermo_params_cond_reaction $log_file $nucleoside_name $component_nucleoside
calc_thermo_params_cond_reaction $RU.log $RU $component_RU
calc_thermo_params_cond_reaction "$TC"_"$ring_conf_fura".log $TC_name $component_TC 

output_other_info $ring_conf_fura $RU

output_results_for_water $environment


cd ../
done

cd ../
done

elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then

for ring_conf_pyra in "${Ring_Conf_pyranose[@]}"; do

cd $ring_conf_pyra

for RU in "${RU_names[@]}"; do

cd $RU

# Output message of status of script
echo "working on $environment/$TC/$ring_conf_pyra/$RU"

# Define if there is one or two or three opt log files for nucleoside
# Define input file for function that calculates the thermo equilibrium
log_file_prefix="$TC"_"$ring_conf_pyra"_"$RU"

# Define input file for function that calculates the thermo equilibrium
unique_opt_file="$log_file_prefix"_min.log
second_opt_file="$log_file_prefix"_2opt.log
third_opt_file="$log_file_prefix"_3opt.log

if [ -f "$third_opt_file" ]; then

log_file=$third_opt_file

echo "3opt" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/n_gaussian_opt.txt

else

if [ -f "$second_opt_file" ]; then

log_file=$second_opt_file

echo "2opt" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/n_gaussian_opt.txt

else

log_file=$unique_opt_file

echo "1opt" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/n_gaussian_opt.txt

fi

fi


# Run function to extract thermo parameters from log files for nucleoside, RU and TC 
component_TC="TC"
component_RU="RU"
component_nucleoside="nucleoside"
TC_name="$TC"_"$ring_conf_pyra"
nucleoside_name="$TC"_"$ring_conf_pyra"_"$RU"
calc_thermo_params_cond_reaction $log_file $nucleoside_name $component_nucleoside
calc_thermo_params_cond_reaction $RU.log $RU $component_RU
calc_thermo_params_cond_reaction "$TC"_"$ring_conf_pyra".log $TC_name $component_TC 

output_other_info $ring_conf_pyra $RU

output_results_for_water $environment


cd ../
done

cd ../
done

elif [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then

for RU in "${RU_names[@]}"; do

cd $RU

# Output message of status of script
echo "working on $environment/$TC/$RU"

# Define if there is one or two or three opt log files for nucleoside
# Define input file for function that calculates the thermo equilibrium
log_file_prefix="$TC"_"$RU"

# Define input file for function that calculates the thermo equilibrium
unique_opt_file="$log_file_prefix"_min.log
second_opt_file="$log_file_prefix"_2opt.log
third_opt_file="$log_file_prefix"_3opt.log

if [ -f "$third_opt_file" ]; then

log_file=$third_opt_file

echo "3opt" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/n_gaussian_opt.txt

else

if [ -f "$second_opt_file" ]; then

log_file=$second_opt_file

echo "2opt" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/n_gaussian_opt.txt

else

log_file=$unique_opt_file

echo "1opt" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/n_gaussian_opt.txt

fi

fi

# Run function to extract thermo parameters from log files for nucleoside, RU and TC 
component_TC="TC"
component_RU="RU"
component_nucleoside="nucleoside"
TC_name="$TC"
nucleoside_name="$TC"_"$RU"
calc_thermo_params_cond_reaction $log_file $nucleoside_name $component_nucleoside
calc_thermo_params_cond_reaction $RU.log $RU $component_RU
calc_thermo_params_cond_reaction "$TC".log $TC_name $component_TC 

ring_conf="none"
output_other_info $ring_conf $RU

output_results_for_water $environment




cd ../
done

elif [ "$TC" == "peptide" ]; then

for RU_ac in "${RU_ac_names[@]}"; do

cd $RU_ac

# Output message of status of script
echo "working on $environment/$TC/$RU_ac"

# Define if there is one or two or three opt log files for nucleoside
# Define input file for function that calculates the thermo equilibrium
log_file_prefix="$TC"_"$RU_ac"

# Define input file for function that calculates the thermo equilibrium
unique_opt_file="$log_file_prefix"_min.log
second_opt_file="$log_file_prefix"_2opt.log
third_opt_file="$log_file_prefix"_3opt.log

if [ -f "$third_opt_file" ]; then

log_file=$third_opt_file

echo "3opt" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/n_gaussian_opt.txt

else

if [ -f "$second_opt_file" ]; then

log_file=$second_opt_file

echo "2opt" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/n_gaussian_opt.txt

else

log_file=$unique_opt_file

echo "1opt" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/condensation_reaction/n_gaussian_opt.txt

fi

fi

# Run function to extract thermo parameters from log files for nucleoside, RU and TC 
component_TC="TC"
component_RU="RU"
component_nucleoside="nucleoside"
TC_name="$TC"
nucleoside_name="$TC"_"$RU_ac"
calc_thermo_params_cond_reaction $log_file $nucleoside_name $component_nucleoside
calc_thermo_params_cond_reaction $RU_ac.log $RU_ac $component_RU
calc_thermo_params_cond_reaction "$TC".log $TC_name $component_TC 

ring_conf="none"
output_other_info $ring_conf $RU_ac

output_results_for_water $environment


cd ../
done

fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT

cd ../


###########################################################################################################################################
############################################## WORK ON CONDENSATION REACTION FOLDER #######################################################
###########################################################################################################################################

# Open folder
cd post_proc_summary_results_csv_graphs/condensation_reaction

# Create and run python script to gen summary csv files and freq hist 
cat >create_csv_files_hist_cond_reaction.py <<'END_SCRIPT'

import sys, os
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from numpy import loadtxt
import seaborn as sb
import scipy as scipy
from scipy import stats
from itertools import repeat
from matplotlib.gridspec import SubplotSpec
import shutil

####################################### Generate necesary variables #############################
#Environment
environment = ['vacuum','water']

# TC
TC_names = ["2deoxy_ribofuranose", "ribofuranose", "threose", "2deoxy_ribopyranose", "ribopyranose", 
            "glycerol", "glyceric_acid", "peptide"]

# RUs
RUs_names=["adenine", "guanine", "cytosine", "thymine", "uracil", "TARC_cbond", "TARC_nbond", "BA_cbond", "BA_nbond", "CA", "melamine"]
RUs_ac_names=["adenine_ac", "guanine_ac", "cytosine_ac", "thymine_ac", "uracil_ac", "TARC_cbond_ac", "TARC_nbond_ac", "BA_cbond_ac", "BA_nbond_ac", "CA_ac", "melamine_ac"]
new_RUs_names = ["A", "G", "C", "T", "U", "TAP-$C^5$", "TAP-N", "BA-$C^5$", "BA-N", "CA", "MM"]
new_RUs_ac_names = ["A", "G", "C", "T", "U", "TAP-$C^5$", "TAP-N", "BA-$C^5$", "BA-N", "CA", "MM"]

# Type RU
type_RU = ["canonical", "canonical", "canonical"]

# Ring conf
Ring_Conf_furanose = ["2endo_alpha", "2endo_beta", "3endo_alpha", "3endo_beta"]
Ring_Conf_pyranose = ["D1C4_alpha", "D1C4_beta", "D4C1_alpha", "D4C1_beta"]

# Thermo param
#thermo_params = ['ΔzpE', 'ΔcorrE', 'ΔH', 'ΔG', 'Δ(TdS)']
thermo_params = ['ΔE', 'ΔE(ZPE)', 'ΔH°', 'ΔG°', 'Δ(TΔS°)']

####################################### Generate necesary functions #############################

##################################################
##                                              ##
## FUNCTIONS TO GENERATE SUBTITLES.             ##
##                                              ##
## Description:                                 ##
## This function generates subtitles.           ##                
##                                              ##
##################################################

def create_subtitle(fig: plt.Figure, grid: SubplotSpec, title: str):
    "Sign sets of subplots with title"
    row = fig.add_subplot(grid)
    # the '\n' is important
    #if title == 'HAsO$_3^-$-2dRibf-RU in vacuum':
    row.set_title(f'{title}\n', fontweight='semibold', fontsize=26)
    #else:
     #   row.set_title(f'{title}', fontweight='semibold', fontsize=26)
    # hide subplot
    row.set_frame_on(False)
    row.axis('off')
    

##################################################
##                                              ##
## FUNCTIONS TO GENERATE FREQ HIST.             ##
##                                              ##
## Description:                                 ##
## This function generates for each thermo      ##
## param and each TC a bar plot of the delta    ##
## param of nucleoside condensation reaction.   ##                   
##                                              ##
##################################################
  
def gen_bar_graph (dataframe, thermo, TC):
    
    # Define number of rows and cols
    rows = 2
    cols = 2
    
    # Create the fig and axn
    if TC == 'peptide':
        fig, axn = plt.subplots(rows, cols, sharex=False, sharey=False, figsize=(24,14))
    else:
        fig, axn = plt.subplots(rows, cols, sharex=False, sharey=False, figsize=(22,12))
        
    # Define the grid
    grid = plt.GridSpec(rows, cols)    
    
    # Add subtitles
    if TC == "2deoxy_ribofuranose":
        create_subtitle(fig, grid[0, ::], '2dRibf-RU in vacuum')
        create_subtitle(fig, grid[1, ::], '2dRibf-RU in water')
    elif TC == "ribofuranose":
        create_subtitle(fig, grid[0, ::], 'Ribf-RU in vacuum')
        create_subtitle(fig, grid[1, ::], 'Ribf-RU in water')
    elif TC == "threose":
        create_subtitle(fig, grid[0, ::], 'Tho-RU in vacuum')
        create_subtitle(fig, grid[1, ::], 'Tho-RU in water')
    elif TC == "2deoxy_ribopyranose":
        create_subtitle(fig, grid[0, ::], '2dRib-RU in vacuum')
        create_subtitle(fig, grid[1, ::], '2dRib-RU in water')
    elif TC == "ribopyranose":
        create_subtitle(fig, grid[0, ::], 'Rib-RU in vacuum')
        create_subtitle(fig, grid[1, ::], 'Rib-RU in water')
    elif TC == "glycerol":
        create_subtitle(fig, grid[0, ::], 'GLY-RU in vacuum')
        create_subtitle(fig, grid[1, ::], 'GLY-RU in water')
    elif TC == "glyceric_acid":
        create_subtitle(fig, grid[0, ::], 'GA-RU in vacuum')
        create_subtitle(fig, grid[1, ::], 'GA-RU in water')
    elif TC == "peptide":
        create_subtitle(fig, grid[0, ::], 'AEG-RU in vacuum')
        create_subtitle(fig, grid[1, ::], 'AEG-RU in water')

    for i, (ax, env, type_ru) in enumerate([(axn.flat[0], "vacuum", "canonical"), 
                                            (axn.flat[1], "vacuum", "non_canonical"), 
                                            (axn.flat[2], "water", "canonical"),
                                            (axn.flat[3], "water", "non_canonical")]):
        
        # Define list of canonical and non canonical RUs
        RUs_canonical = ["A", "G", "C", "T", "U"]
        RUs_non_canonical = ["TAP-$C^5$", "TAP-N", "BA-$C^5$", "BA-N", "CA", "MM"]
        RUs_canonical_ac = ["A", "G", "C", "T", "U"]
        RUs_non_canonical_ac = ["TAP-$C^5$", "TAP-N", "BA-$C^5$", "BA-N", "CA", "MM"]
        
        
        # Obtain rows from database that belong to vacuum or water env
        env_df = dataframe[dataframe['Env'].isin([env])]
        
        # From this new df obtain rows for either canonical or non canonical bases
        if type_ru == "canonical" and TC != "peptide":
            ru_df = env_df[env_df['RUs'].isin(RUs_canonical)]
        elif type_ru == "non_canonical" and TC != "peptide":
            ru_df = env_df[env_df['RUs'].isin(RUs_non_canonical)]
        elif type_ru == "canonical" and TC == "peptide":
            ru_df = env_df[env_df['RUs'].isin(RUs_canonical_ac)]
        elif type_ru == "non_canonical" and TC == "peptide":
            ru_df = env_df[env_df['RUs'].isin(RUs_non_canonical_ac)]
        
        # Set indexes of database
        ru_df = ru_df.set_index("RUs")
        #print(ru_df)
    
        # Generate bar graph now
        if env == "vacuum" and type_ru == "canonical":
            if TC in ['2deoxy_ribofuranose', 'ribofuranose', 'threose', '2deoxy_ribopyranose', 'ribopyranose']:
                ax1 = ru_df.plot(kind='bar', ax=ax, colormap='Greys', rot=360, width=0.8, edgecolor = "black", 
                                 xlabel='', ylabel='ΔG (kJ/mol)', legend=False, fontsize=20)
            else:
                ax1 = ru_df.plot(kind='bar', ax=ax, color='Grey', rot=360, width=0.8, edgecolor = "black", 
                                 xlabel='', ylabel='ΔG (kJ/mol)', legend=False, fontsize=20)
            ax1.set_xticks([]) 
            ax1.set_ylabel(str(thermo)+ " (kJ/mol)", fontsize=20, fontdict=dict(weight='bold'))
            ax1.set_title('Canonical RUs', weight='bold', fontsize=22)
            for container in ax1.containers:
                ax1.bar_label(container, fmt='%.1f', fontsize=14)
        elif env == "vacuum" and type_ru == "non_canonical":
            if TC in ['2deoxy_ribofuranose', 'ribofuranose', 'threose', '2deoxy_ribopyranose', 'ribopyranose']:
                ax2 = ru_df.plot(kind='bar', ax=ax, colormap='Greys', rot=360, width=0.8, edgecolor = "black", 
                                 xlabel='', legend=False, fontsize=20)
            else:
                ax2 = ru_df.plot(kind='bar', ax=ax, color='Grey', rot=360, width=0.8, edgecolor = "black", 
                                 xlabel='', legend=False, fontsize=20)
            ax2.set_xticks([])
            ax2.set_title('Non-canonical RUs', weight='bold', fontsize=22)
            for container in ax2.containers:
                ax2.bar_label(container, fmt='%.1f', fontsize=14)
        elif env == "water" and type_ru == "canonical":
            if TC in ['2deoxy_ribofuranose', 'ribofuranose', 'threose', '2deoxy_ribopyranose', 'ribopyranose']:
                ax3 = ru_df.plot(kind='bar', ax=ax, colormap='Greys', rot=360, width=0.8, edgecolor = "black", 
                                 xlabel='', ylabel='ΔG (kJ/mol)', legend=False, fontsize=20)
            else:
                ax3 = ru_df.plot(kind='bar', ax=ax, color='Grey', rot=360, width=0.8, edgecolor = "black", 
                                 xlabel='', ylabel='ΔG (kJ/mol)', legend=False, fontsize=20)
            ax3.set_ylabel(str(thermo) + " (kJ/mol)", fontsize=20, fontdict=dict(weight='bold'))
            #ax3.set_title('Water canonical RUs', weight='bold', fontsize=20)
            for container in ax3.containers:
                ax3.bar_label(container, fmt='%.1f', fontsize=14)
        elif env == "water" and type_ru == "non_canonical":
            if TC in ['2deoxy_ribofuranose', 'ribofuranose', 'threose', '2deoxy_ribopyranose', 'ribopyranose']:
                ax4 = ru_df.plot(kind='bar', ax=ax, colormap='Greys', rot=360, width=0.8, edgecolor = "black", 
                                 xlabel='', legend=True, fontsize=20)
            else:
                ax4 = ru_df.plot(kind='bar', ax=ax, color='Grey', rot=360, width=0.8, edgecolor = "black", 
                                 xlabel='', fontsize=20, legend=False)
                #ax4.legend([str(thermo)+"_react"])
            #ax4.set_title('Water non_canonical RUs', weight='bold', fontsize=22)
            if TC in ['2deoxy_ribofuranose', 'ribofuranose', 'threose']:
                ax4.legend(['α-South', 'β-South', 'α-North', 'β-North'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=16)
            elif TC in ['2deoxy_ribopyranose', 'ribopyranose']:
                ax4.legend(['α-$^1C_4$', 'β-$^1C_4$', 'α-$^4C_1$', 'β-$^4C_1$'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=16)
            #else:
                #ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=16)
            for container in ax4.containers:
                ax4.bar_label(container, fmt='%.1f', fontsize=14)
            
            
    # Set parameters for axis labels and titles       
    for axis in [ax1, ax2, ax3, ax4]:
        axis.axhline(y = 0.0, color = 'black', linestyle = '--')
        labels = axis.get_xticklabels() 
        #for label in labels:
            #label.set_fontweight('bold')
    
    plt.tight_layout()
    plt.savefig('Condensation_react_'+str(TC)+'_'+str(thermo)+'.pdf')
    
 
    

####################################### Create general database #################################

# Load files with different info in variables
RU = np.loadtxt("RUs_names_list.txt", dtype='str')
TC = np.loadtxt("TCs_names_list.txt", dtype='str')
initial_ring_conf = np.loadtxt("ring_conf_names_list.txt", dtype='str')
environment = np.loadtxt("environment_list.txt", dtype='str')

# Load energies txt files for each component in different variable
zpE_water = np.loadtxt("water_zero_point_energies.txt", dtype='float')
corrE_water = np.loadtxt("water_corrected_energies.txt", dtype='float')
H_water = np.loadtxt("water_enthalpies.txt", dtype='float')
G_water = np.loadtxt("water_free_energies.txt", dtype='float')

zpE_nucleoside = np.loadtxt("nucleosides_zero_point_energies.txt", dtype='float')
corrE_nucleoside = np.loadtxt("nucleosides_corrected_energies.txt", dtype='float')
H_nucleoside = np.loadtxt("nucleosides_enthalpies.txt", dtype='float')
G_nucleoside = np.loadtxt("nucleosides_free_energies.txt", dtype='float')

zpE_RU = np.loadtxt("RUs_zero_point_energies.txt", dtype='float')
corrE_RU = np.loadtxt("RUs_corrected_energies.txt", dtype='float')
H_RU = np.loadtxt("RUs_enthalpies.txt", dtype='float')
G_RU = np.loadtxt("RUs_free_energies.txt", dtype='float')

zpE_TC = np.loadtxt("TCs_zero_point_energies.txt", dtype='float')
corrE_TC = np.loadtxt("TCs_corrected_energies.txt", dtype='float')
H_TC = np.loadtxt("TCs_enthalpies.txt", dtype='float')
G_TC = np.loadtxt("TCs_free_energies.txt", dtype='float')

db_nucleosides_cond_reaction = pd.DataFrame({'Environment': environment, 'TC': TC, 'Initial ring conf': initial_ring_conf, 'RU': RU, 
     'Initial_sugar_ring_conf': initial_ring_conf, 'zpE_nucleoside(kJ/mol)': zpE_nucleoside, 
     'corrE_nucleoside(kJ/mol)': corrE_nucleoside, 'H_nucleoside(kJ/mol)': H_nucleoside, 
     'G_nucleoside(kJ/mol)': G_nucleoside, 'zpE_RU(kJ/mol)': zpE_RU, 
     'corrE_RU(kJ/mol)': corrE_RU, 'H_RU(kJ/mol)': H_RU, 
     'G_RU(kJ/mol)': G_RU, 'zpE_TC(kJ/mol)': zpE_TC, 'corrE_TC(kJ/mol)': corrE_TC, 'H_TC(kJ/mol)': H_TC,   
     'G_TC(kJ/mol)': G_TC, 'zpE_water(kJ/mol)': zpE_water, 'corrE_water(kJ/mol)': corrE_water,
     'H_water(kJ/mol)': H_water, 'G_water(kJ/mol)': G_water})
     
# Obtain the TdS for the different components   
db_nucleosides_cond_reaction['TdS_nucleoside(kJ/mol)'] = db_nucleosides_cond_reaction['H_nucleoside(kJ/mol)'] - db_nucleosides_cond_reaction['G_nucleoside(kJ/mol)']
db_nucleosides_cond_reaction['TdS_RU(kJ/mol)'] = db_nucleosides_cond_reaction['H_RU(kJ/mol)'] - db_nucleosides_cond_reaction['G_RU(kJ/mol)']
db_nucleosides_cond_reaction['TdS_TC(kJ/mol)'] = db_nucleosides_cond_reaction['H_TC(kJ/mol)'] - db_nucleosides_cond_reaction['G_TC(kJ/mol)']
db_nucleosides_cond_reaction['TdS_water(kJ/mol)'] = db_nucleosides_cond_reaction['H_water(kJ/mol)'] - db_nucleosides_cond_reaction['G_water(kJ/mol)']

# Obtain the Δ thermo parameters 
db_nucleosides_cond_reaction['ΔE_reaction(kJ/mol)'] = (db_nucleosides_cond_reaction['zpE_water(kJ/mol)'] + db_nucleosides_cond_reaction['zpE_nucleoside(kJ/mol)']) - (db_nucleosides_cond_reaction['zpE_RU(kJ/mol)'] + db_nucleosides_cond_reaction['zpE_TC(kJ/mol)'])

db_nucleosides_cond_reaction['ΔE(ZPE)_reaction(kJ/mol)'] = (db_nucleosides_cond_reaction['corrE_water(kJ/mol)'] + db_nucleosides_cond_reaction['corrE_nucleoside(kJ/mol)']) - (db_nucleosides_cond_reaction['corrE_RU(kJ/mol)'] + db_nucleosides_cond_reaction['corrE_TC(kJ/mol)'])

db_nucleosides_cond_reaction['ΔH°_reaction(kJ/mol)'] = (db_nucleosides_cond_reaction['H_water(kJ/mol)'] + db_nucleosides_cond_reaction['H_nucleoside(kJ/mol)']) - (db_nucleosides_cond_reaction['H_RU(kJ/mol)'] + db_nucleosides_cond_reaction['H_TC(kJ/mol)'])

db_nucleosides_cond_reaction['ΔG°_reaction(kJ/mol)'] = (db_nucleosides_cond_reaction['G_water(kJ/mol)'] + db_nucleosides_cond_reaction['G_nucleoside(kJ/mol)']) - (db_nucleosides_cond_reaction['G_RU(kJ/mol)'] + db_nucleosides_cond_reaction['G_TC(kJ/mol)'])

db_nucleosides_cond_reaction['Δ(TΔS°)_reaction(kJ/mol)'] = (db_nucleosides_cond_reaction['TdS_water(kJ/mol)'] + db_nucleosides_cond_reaction['TdS_nucleoside(kJ/mol)']) - (db_nucleosides_cond_reaction['TdS_RU(kJ/mol)'] + db_nucleosides_cond_reaction['TdS_TC(kJ/mol)'])

# Convert pandas dataframe to csv file
db_nucleosides_cond_reaction.to_csv('condensation_reaction_allData_nucleosides.csv', index=False)

####################################################################################################################

# Read the csv file with all info
db = pd.read_csv ('condensation_reaction_allData_nucleosides.csv')

# Change name of RUs in the field of this final dataframe
for i, RU in enumerate(RUs_names):
    #db.replace({'RU': RU}, {'RU': new_RUs_names[i]}, regex=True)
    # Change name of RUs in the field of this final dataframe
    db['RU'] = db['RU'].str.replace(RU,new_RUs_names[i])
    
for i, RU_ac in enumerate(RUs_ac_names):
    # Change name of RUs in the field of this final dataframe
    db['RU'] = db['RU'].str.replace(RU_ac,new_RUs_ac_names[i])
    

####################################### Loop through different thermo and TCs and generate graphs and csv files #################################
a = ['vacuum','water']
n = sorted (a*11)
	
# Select columns that include only thermo param
for thermo in thermo_params:
    db_to_work = db[["Environment", "TC", "Initial ring conf", "RU", str(thermo)+"_reaction(kJ/mol)"]]
    
    # Select from previous db columns for specific TC
    for TC in TC_names:
    
        TC_df = db_to_work[db_to_work['TC'] == TC]
        
        # Work on furanose sugars
        if TC in ['2deoxy_ribofuranose', 'ribofuranose', 'threose']: # FURANOSES
            
            # Generate final database TC, RU and env fields
            env_list = n
            #env_list = sorted (env_list)
            RUs_list = new_RUs_names * 2
            final_db = pd.DataFrame({'Env': env_list, 'TCs': TC, 'RUs': RUs_list})
            
            # Output entries for specific ring conformation
            for ring_conf in Ring_Conf_furanose: 
                db_ring_conf = TC_df[TC_df['Initial ring conf'].str.contains(ring_conf)]

                # Add previous field to final db
                final_db[str(ring_conf)+'_'+str(thermo)+'_react'] = db_ring_conf[str(thermo)+'_reaction(kJ/mol)'].to_numpy()
                       
            # Obtain freq hist
            gen_bar_graph (final_db, thermo, TC)
            final_db.to_csv('summary_'+str(thermo)+'_'+str(TC)+'.csv', index=False)
        
        elif TC in ["2deoxy_ribopyranose", "ribopyranose"]: ## PYRANOSES
            
            # Generate final database TC, RU and env fields
            env_list = n
            #env_list = sorted (environment*11)
            RUs_list = new_RUs_names * 2
            final_db = pd.DataFrame({'Env': env_list, 'TCs': TC, 'RUs': RUs_list})
            
            # Output entries for specific ring conformation
            for ring_conf in Ring_Conf_pyranose: 
                db_ring_conf = TC_df[TC_df['Initial ring conf'].str.contains(ring_conf)]

                # Add previous field to final db
                final_db[str(ring_conf)+'_'+str(thermo)+'_react'] = db_ring_conf[str(thermo)+'_reaction(kJ/mol)'].to_numpy()
            
            # Obtain freq hist
            gen_bar_graph (final_db, thermo, TC)
            final_db.to_csv('summary_'+str(thermo)+'_'+str(TC)+'.csv', index=False)
            
        elif TC in ["glycerol", "glyceric_acid"]: ## NON SUGAR
            
            # Generate final database TC, RU and env fields
            env_list = n
            #env_list = sorted (environment*11)
            RUs_list = new_RUs_names * 2
            final_db = pd.DataFrame({'Env': env_list, 'TCs': TC, 'RUs': RUs_list})
        
            db_ring_conf = TC_df[TC_df['Initial ring conf'] == 'none']
            
            # Add previous field to final db
            final_db[str(thermo)+'_react'] = db_ring_conf[str(thermo)+'_reaction(kJ/mol)'].to_numpy()
            
            # Obtain freq hist     
            gen_bar_graph (final_db, thermo, TC)
            final_db.to_csv('summary_'+str(thermo)+'_'+str(TC)+'.csv', index=False)
            
        else: ## PEPTIDE
        
            # Generate final database TC, RU and env fields
            #env_list = sorted(environment * 11)
            env_list = n
            #env_list = sorted (env_list)
            RUs_list = new_RUs_ac_names * 2
            final_db = pd.DataFrame({'Env': env_list, 'TCs': TC, 'RUs': RUs_list})
        
            db_ring_conf = TC_df[TC_df['Initial ring conf'] == 'none']
                
            # Add previous field to final db
            final_db[str(thermo)+'_react'] = db_ring_conf[str(thermo)+'_reaction(kJ/mol)'].to_numpy()
            
            # Obtain freq hist     
            gen_bar_graph (final_db, thermo, TC)
            final_db.to_csv('summary_'+str(thermo)+'_'+str(TC)+'.csv', index=False)

# Move files to their new folders
for TC in TC_names:
    os.mkdir(TC)
    for thermo in thermo_params:
        shutil.move('Condensation_react_'+str(TC)+'_'+str(thermo)+'.pdf', str(TC)+'/Condensation_react_'+str(TC)+'_'+str(thermo)+'.pdf')
        shutil.move('summary_'+str(thermo)+'_'+str(TC)+'.csv', str(TC)+'/summary_'+str(thermo)+'_'+str(TC)+'.csv') 
                    
END_SCRIPT

# Run script.py with python3
python3 create_csv_files_hist_cond_reaction.py  


















