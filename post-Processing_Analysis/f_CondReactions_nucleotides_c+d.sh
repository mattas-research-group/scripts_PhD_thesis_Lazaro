#!/bin/bash


#	f_CondReactions_nucleotides_c+d.sh

#	18-jun-2023	Edited options for subtitle and group files in TC folders in python function
# 	04-mar-2023	Created

# Description:
# This script obtains the delta thermo params for the 
# condensation reaction of process a+b for nucleotides.

# TO RUN
# ./f_CondReactions_nucleotides_c+d.sh  

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a RU_names=("adenine" "guanine" "cytosine" "thymine" "uracil" "TARC_cbond" "TARC_nbond" "BA_cbond" "BA_nbond" "CA" "melamine")
declare -a IL_names=("phosphate" "arsenate")
declare -a thermo_names=("zeropoint_energy_kj_mol" "corrected_energy_kj_mol" "G_kj_mol" "entalphy_kj_mol")
declare -a thermo_names_ru=("zeropoint_energy_kcal_mol" "corrected_energy_kcal_mol" "G_kcal_mol" "entalphy_kcal_mol")
declare -a thermos=("zero_point_energies" "corrected_energies" "free_energies" "enthalpies")




function num_gaussian_opt {

# Assign variables names to entries of function
local log_file_prefix=$1

# Run if then statement to define maximum number of opt+freq
find -type f -name "$log_file_prefix*" -and -name "*min.log" > list_gaussian_min.txt
find -type f -name "$log_file_prefix*" -and -name "*opt.log" > list_gaussian_opt.txt
sort -n list_gaussian_opt.txt > tmp.txt && mv tmp.txt list_gaussian_opt.txt
cat list_gaussian_min.txt list_gaussian_opt.txt > final_list_gaussian.txt

# Delete first two characters ./ from each line
cut -c 3- final_list_gaussian.txt > tmp.txt && mv tmp.txt final_list_gaussian.txt

# Select last file to work on it
log_file=$( tail -n 1 final_list_gaussian.txt )

echo $log_file > log_file_from_nucleotide_working_on.txt

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
local IL=$4

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
if [ "$component_type" == "nucleotide" ]; then
echo $zero_point_energy >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/nucleotides_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/nucleotides_corrected_energies.txt
echo $enthalpy >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/nucleotides_enthalpies.txt
echo $free_energy >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/nucleotides_free_energies.txt

fi

}







function output_other_info {

# Declare variables
local environment=$1
local TC=$2
local ring_conf=$3
local IL=$4
local RU=$5

# Output names of components in general files
echo $environment >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/environment_list.txt
echo $TC >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/TCs_names_list.txt
echo $ring_conf >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/ring_conf_names_list.txt
echo $IL >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/ILs_names_list.txt
echo $RU >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/RUs_names_list.txt

# Output info on thermo params of TC, IL, backbone and RU
for i in "${!thermo_names[@]}"; do      

# Obtain info for IL
IL_thermo=`cat "$IL"_"${thermo_names[$i]}".txt`
echo $IL_thermo >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/ILs_"${thermos[$i]}".txt

# Obtain info for TC and backbone
if [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then

# Output info for TC
TC_thermo=`cat "$TC"_"${thermo_names[$i]}".txt`
echo $TC_thermo >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/TCs_"${thermos[$i]}".txt

# Output info for backbone
backbone_thermo=`cat "$TC"_"$IL"_"${thermo_names[$i]}".txt`
echo $backbone_thermo >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/backbone_"${thermos[$i]}".txt
#echo "$TC"_"$RU"_"$IL"_"${thermo_names[$i]}".txt >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/nucleotides_"${thermos[$i]}".txt

else

# Output info for TC
TC_thermo=`cat "$TC"_"$ring_conf"_"${thermo_names[$i]}".txt`
echo $TC_thermo >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/TCs_"${thermos[$i]}".txt

# Output info for backbone
backbone_thermo=`cat "$TC"_"$ring_conf"_"$IL"_"${thermo_names[$i]}".txt`
echo $backbone_thermo >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/backbone_"${thermos[$i]}".txt
#echo "$TC"_"$ring_conf"_"$RU"_"$IL"_"${thermo_names[$i]}".txt >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/nucleotides_"${thermos[$i]}".txt

fi

# Output info for the RU
RU_thermo=`cat "$RU"_"${thermo_names_ru[$i]}".txt`
echo $RU_thermo >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/RUs_"${thermos[$i]}".txt

done

}







function output_results_for_water {

local environment=$1

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

echo $zpE_water >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/water_zero_point_energies.txt
echo $corrE_water >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/water_corrected_energies.txt
echo $H_water >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/water_enthalpies.txt
echo $G_water >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL/water_free_energies.txt

}






##########################################################################################################################
############################################## WORK ON EACH FOLDER #######################################################
##########################################################################################################################

# Open corresponding folder
mkdir nucleotides_c+d/post_proc_summary_results_csv_graphs
mkdir nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction

# Create IL folder inside the post-processing folder
for IL in "${IL_names[@]}"; do
mkdir nucleotides_c+d/post_proc_summary_results_csv_graphs/condensation_reaction/$IL
done




cd nucleotides_c+d/final_gaussian_opt

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

# Define input file for function that calculates the number of opt+freq log files
log_file_prefix="$TC"_"$ring_conf_fura"_"$RU"_"$IL"

# Run function to determine how many opt+freq were runned and get log file
num_gaussian_opt $log_file_prefix

# Run function that extract all the energies from this file and print them in corresponding summary files in processing folder
calc_thermo_params_cond_reaction $log_file $log_file_prefix "nucleotide" $IL


# Output name of components
output_other_info $environment $TC $ring_conf_fura $IL $RU

# Output energies from water
output_results_for_water $environment





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

# Define input file for function that calculates the number of opt+freq log files
log_file_prefix="$TC"_"$ring_conf_pyra"_"$RU"_"$IL"

# Run function to determine how many opt+freq were runned and get log file
num_gaussian_opt $log_file_prefix

# Run function that extract all the energies from this file and print them in corresponding summary files in processing folder
calc_thermo_params_cond_reaction $log_file $log_file_prefix "nucleotide" $IL


# Output name of components
output_other_info $environment $TC $ring_conf_pyra $IL $RU

# Output energies from water
output_results_for_water $environment







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

# Define input file for function that calculates the number of opt+freq log files
log_file_prefix="$TC"_"$RU"_"$IL"

# Run function to determine how many opt+freq were runned and get log file
num_gaussian_opt $log_file_prefix

# Run function that extract all the energies from this file and print them in corresponding summary files in processing folder
calc_thermo_params_cond_reaction $log_file $log_file_prefix "nucleotide" $IL


# Output name of components
output_other_info $environment $TC $ring_conf $IL $RU

# Output energies from water
output_results_for_water $environment







cd ../
done

cd ../
done

fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT

cd ../ # LEAVE FINAL_GAUSSIAN_OPT FOLDER

###########################################################################################################################################
############################################## WORK ON CONDENSATION REACTION FOLDER #######################################################
###########################################################################################################################################
pdw 

# Open folder
cd post_proc_summary_results_csv_graphs/condensation_reaction

for IL in "${IL_names[@]}"; do

cd $IL

export IL

# Create and run python script to gen summary csv files and freq hist 
cat >create_csv_files_hist_cond_reaction.py <<'END_SCRIPT'

################## Import necesary libraries ##################
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

################ Create necesary stacked lists #################
#Environment
environment = ["vacuum", "water"]

# Import environmental variables from bash script for angles in degrees
IL_name = os.environ.get("IL")

# TC
TC_names = ["2deoxy_ribofuranose", "ribofuranose", "threose", "2deoxy_ribopyranose", "ribopyranose", 
            "glycerol", "glyceric_acid"]

# RUs
RUs_names=["adenine", "guanine", "cytosine", "thymine", "uracil", "TARC_cbond", "TARC_nbond", "BA_cbond", "BA_nbond", "CA", "melamine"]
new_RUs_names = ["A", "G", "C", "T", "U", "TAP-$C^5$", "TAP-N", "BA-$C^5$", "BA-N", "CA", "MM"]

# Type RU
type_RU = ["canonical", "canonical", "canonical"]

# Ring conf
Ring_Conf_furanose = ["2endo_alpha", "2endo_beta", "3endo_alpha", "3endo_beta"]
new_ring_conf_fura = ['α-$^2T_3$', 'β-$^2T_3$', 'α-$^3T_2$', 'β-$^3T_2$']
Ring_Conf_pyranose = ["D1C4_alpha", "D1C4_beta", "D4C1_alpha", "D4C1_beta"]
new_ring_conf_pyra = ['α-$^1C_4$', 'β-$^1C_4$', 'α-$^4C_1$', 'β-$^4C_1$']

# Thermo param
thermo_params = ['ΔE', 'ΔE(ZPE)', 'ΔH°', 'ΔG°', 'Δ(TΔS°)']

# Generate necesary functions

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
## param of reaction a+b.                       ##                   
##                                              ##
##################################################

def gen_bar_graph (dataframe, thermo, TC, IL_name):
    
    fig, axn = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(22,12)) 

    for i, (ax, env, type_ru) in enumerate([(axn.flat[0], "vacuum", "canonical"), 
                                            (axn.flat[1], "vacuum", "non_canonical"), 
                                            (axn.flat[2], "water", "canonical"),
                                            (axn.flat[3], "water", "non_canonical")]):
        
        # Define list of canonical and non canonical RUs
        RUs_canonical = ["A", "G", "C", "T", "U"]
        RUs_non_canonical = ["TAP-$C^5$", "TAP-N", "BA-$C^5$", "BA-N", "CA", "MM"]
        
        # Obtain rows from database that belong to vacuum or water env
        env_df = dataframe[dataframe['Env'].isin([env])]
        
        # From this new df obtain rows for either canonical or non canonical bases
        if type_ru == "canonical" and TC != "peptide":
            ru_df = env_df[env_df['RUs'].isin(RUs_canonical)]
        elif type_ru == "non_canonical" and TC != "peptide":
            ru_df = env_df[env_df['RUs'].isin(RUs_non_canonical)]
              
        # Set indexes of database
        ru_df = ru_df.set_index("RUs")
    
        # Generate bar graph now
        if env == "vacuum" and type_ru == "canonical":
            if TC in ['2deoxy_ribofuranose', 'ribofuranose', 'threose', '2deoxy_ribopyranose', 'ribopyranose']:
                ax1 = ru_df.plot(kind='bar', ax=ax, colormap='Greys', rot=360, width=0.6, edgecolor = "black", 
                                 xlabel='', ylabel='ΔG (kcal/mol)', legend=False, fontsize=28)
            else:
                ax1 = ru_df.plot(kind='bar', ax=ax, color='Grey', rot=360, width=0.6, edgecolor = "black", 
                                 xlabel='', ylabel='ΔG (kcal/mol)', legend=False, fontsize=28)
            ax1.set_xticks([]) 
            ax1.set_ylabel(str(thermo)+"$_{c+d}$ (kJ/mol)", fontsize=28, fontdict=dict(weight='bold'))
            ax1.set_title('Canonical RUs', weight='bold', fontsize=32)
            for container in ax1.containers:
                ax1.bar_label(container, fmt='%.1f', fontsize=18)
            # Generate vertical lines
            [ax1.axvline(x, color = 'grey', linestyle='--') for x in [0.5, 1.5, 2.5, 3.5]]
        elif env == "vacuum" and type_ru == "non_canonical":
            if TC in ['2deoxy_ribofuranose', 'ribofuranose', 'threose', '2deoxy_ribopyranose', 'ribopyranose']:
                ax2 = ru_df.plot(kind='bar', ax=ax, colormap='Greys', rot=360, width=0.6, edgecolor = "black", 
                                 xlabel='', legend=False, fontsize=28)
            else:
                ax2 = ru_df.plot(kind='bar', ax=ax, color='Grey', rot=360, width=0.6, edgecolor = "black", 
                                 xlabel='', legend=False, fontsize=28)
            ax2.set_xticks([])
            ax2.set_title('Non-canonical RUs', weight='bold', fontsize=32)
            for container in ax2.containers:
                ax2.bar_label(container, fmt='%.1f', fontsize=18)
            # Generate vertical lines
            [ax2.axvline(x, color = 'grey', linestyle='--') for x in [0.5, 1.5, 2.5, 3.5, 4.5]] 
        elif env == "water" and type_ru == "canonical":
            if TC in ['2deoxy_ribofuranose', 'ribofuranose', 'threose', '2deoxy_ribopyranose', 'ribopyranose']:
                ax3 = ru_df.plot(kind='bar', ax=ax, colormap='Greys', rot=360, width=0.6, edgecolor = "black", 
                                 xlabel='', ylabel='ΔG (kJ/mol)', legend=False, fontsize=28)
            else:
                ax3 = ru_df.plot(kind='bar', ax=ax, color='Grey', rot=360, width=0.6, edgecolor = "black", 
                                 xlabel='', ylabel='ΔG (kJ/mol)', legend=False, fontsize=28)
            ax3.set_ylabel(str(thermo)+"$_{c+d}$ (kJ/mol)", fontsize=28, fontdict=dict(weight='bold'))
            #ax3.set_title('Water canonical RUs', weight='bold', fontsize=18)
            for container in ax3.containers:
                ax3.bar_label(container, fmt='%.1f', fontsize=18)
            # Generate vertical lines
            [ax3.axvline(x, color = 'grey', linestyle='--') for x in [0.5, 1.5, 2.5, 3.5]] 
        elif env == "water" and type_ru == "non_canonical":
            if TC in ['2deoxy_ribofuranose', 'ribofuranose', 'threose', '2deoxy_ribopyranose', 'ribopyranose']:
                ax4 = ru_df.plot(kind='bar', ax=ax, colormap='Greys', rot=360, width=0.6, edgecolor = "black", 
                                 xlabel='', legend=False, fontsize=28)
            else:
                ax4 = ru_df.plot(kind='bar', ax=ax, color='Grey', rot=360, width=0.6, edgecolor = "black", 
                                 xlabel='', fontsize=28, legend=False)
                #ax4.legend([str(thermo)+"_react"])
            # Generate vertical lines
            [ax4.axvline(x, color = 'grey', linestyle='--') for x in [0.5, 1.5, 2.5, 3.5, 4.5]]
            # Set legend
            if TC in ['2deoxy_ribofuranose', 'ribofuranose', 'threose', '2deoxy_ribopyranose', 'ribopyranose']:
                ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=24)
            for container in ax4.containers:
                ax4.bar_label(container, fmt='%.1f', fontsize=18)
            
    # Set parameters for axis labels and titles       
    for axis in [ax1, ax2, ax3, ax4]:
        axis.axhline(y = 0.0, color = 'black', linestyle = '--')
        labels = axis.get_xticklabels() 
        
    
    plt.tight_layout()
    plt.savefig('Condensation_react_'+str(TC)+'_'+str(thermo)+'_'+str(IL_name)+'.pdf')
    
####################################### Create general database #################################
# Load files with different info in variables
RU = np.loadtxt("RUs_names_list.txt", dtype='str')
TC = np.loadtxt("TCs_names_list.txt", dtype='str')
initial_ring_conf = np.loadtxt("ring_conf_names_list.txt", dtype='str')
environment = np.loadtxt("environment_list.txt", dtype='str')
IL = np.loadtxt("ILs_names_list.txt", dtype='str')

# Load energies txt files for each component in different variable
# Nucleotide
zpE_nucleotide = np.loadtxt("nucleotides_zero_point_energies.txt", dtype='float')
corrE_nucleotide = np.loadtxt("nucleotides_corrected_energies.txt", dtype='float')
H_nucleotide = np.loadtxt("nucleotides_enthalpies.txt", dtype='float')
G_nucleotide = np.loadtxt("nucleotides_free_energies.txt", dtype='float')

# Backbone
zpE_backbone = np.loadtxt("backbone_zero_point_energies.txt", dtype='float')
corrE_backbone = np.loadtxt("backbone_corrected_energies.txt", dtype='float')
H_backbone = np.loadtxt("backbone_enthalpies.txt", dtype='float')
G_backbone = np.loadtxt("backbone_free_energies.txt", dtype='float')

# IL
zpE_IL = np.loadtxt("ILs_zero_point_energies.txt", dtype='float')
corrE_IL = np.loadtxt("ILs_corrected_energies.txt", dtype='float')
H_IL = np.loadtxt("ILs_enthalpies.txt", dtype='float')
G_IL = np.loadtxt("ILs_free_energies.txt", dtype='float')

# TC
zpE_TC = np.loadtxt("TCs_zero_point_energies.txt", dtype='float')
corrE_TC = np.loadtxt("TCs_corrected_energies.txt", dtype='float')
H_TC = np.loadtxt("TCs_enthalpies.txt", dtype='float')
G_TC = np.loadtxt("TCs_free_energies.txt", dtype='float')

# RU
zpE_RU = np.loadtxt("RUs_zero_point_energies.txt", dtype='float')
corrE_RU = np.loadtxt("RUs_corrected_energies.txt", dtype='float')
H_RU = np.loadtxt("RUs_enthalpies.txt", dtype='float')
G_RU = np.loadtxt("RUs_free_energies.txt", dtype='float')

# Water
zpE_water = np.loadtxt("water_zero_point_energies.txt", dtype='float')
corrE_water = np.loadtxt("water_corrected_energies.txt", dtype='float')
H_water = np.loadtxt("water_enthalpies.txt", dtype='float')
G_water = np.loadtxt("water_free_energies.txt", dtype='float')

# Create pandas dataframe
db_nucleotides_cond_reaction = pd.DataFrame({'Environment': environment, 
                                             'TC': TC,
                                             'RU': RU, 
                                             'Initial_sugar_ring_conf': initial_ring_conf, 
                                             'IL': IL, 
                                             'zpE_nucleotide(kJ/mol)': zpE_nucleotide, 
                                             'corrE_nucleotide(kJ/mol)': corrE_nucleotide, 
                                             'H_nucleotide(kJ/mol)': H_nucleotide, 
                                             'G_nucleotide(kJ/mol)': G_nucleotide,
                                             'zpE_backbone(kJ/mol)': zpE_backbone, 
                                             'corrE_backbone(kJ/mol)': corrE_backbone, 
                                             'H_backbone(kJ/mol)': H_backbone, 
                                             'G_backbone(kJ/mol)': G_backbone, 
                                             'zpE_IL(kJ/mol)': zpE_IL, 
                                             'corrE_IL(kJ/mol)': corrE_IL, 
                                             'H_IL(kJ/mol)': H_IL, 
                                             'G_IL(kJ/mol)': G_IL, 
                                             'zpE_RU(kJ/mol)': zpE_RU, 
                                             'corrE_RU(kJ/mol)': corrE_RU, 
                                             'H_RU(kJ/mol)': H_RU, 
                                             'G_RU(kJ/mol)': G_RU, 
                                             'zpE_TC(kJ/mol)': zpE_TC, 
                                             'corrE_TC(kJ/mol)': corrE_TC, 
                                             'H_TC(kJ/mol)': H_TC, 
                                             'G_TC(kJ/mol)': G_TC, 
                                             'zpE_water(kJ/mol)': zpE_water, 
                                             'corrE_water(kJ/mol)': corrE_water, 
                                             'H_water(kJ/mol)': H_water, 
                                             'G_water(kJ/mol)': G_water})
                                             
# Obtain the TdS for the different components   
db_nucleotides_cond_reaction['TdS_nucleotide(kJ/mol)'] = db_nucleotides_cond_reaction['H_nucleotide(kJ/mol)'] - db_nucleotides_cond_reaction['G_nucleotide(kJ/mol)']
db_nucleotides_cond_reaction['TdS_backbone(kJ/mol)'] = db_nucleotides_cond_reaction['H_backbone(kJ/mol)'] - db_nucleotides_cond_reaction['G_backbone(kJ/mol)']
db_nucleotides_cond_reaction['TdS_IL(kJ/mol)'] = db_nucleotides_cond_reaction['H_IL(kJ/mol)'] - db_nucleotides_cond_reaction['G_IL(kJ/mol)']
db_nucleotides_cond_reaction['TdS_RU(kJ/mol)'] = db_nucleotides_cond_reaction['H_RU(kJ/mol)'] - db_nucleotides_cond_reaction['G_RU(kJ/mol)']
db_nucleotides_cond_reaction['TdS_TC(kJ/mol)'] = db_nucleotides_cond_reaction['H_TC(kJ/mol)'] - db_nucleotides_cond_reaction['G_TC(kJ/mol)']
db_nucleotides_cond_reaction['TdS_water(kJ/mol)'] = db_nucleotides_cond_reaction['H_water(kJ/mol)'] - db_nucleotides_cond_reaction['G_water(kJ/mol)']

# Obtain the Δ thermo parameters 
# For dzpE
db_nucleotides_cond_reaction['ΔzpE_reaction_c(kJ/mol)'] = (db_nucleotides_cond_reaction['zpE_water(kJ/mol)'] + db_nucleotides_cond_reaction['zpE_backbone(kJ/mol)']) - (db_nucleotides_cond_reaction['zpE_IL(kJ/mol)'] + db_nucleotides_cond_reaction['zpE_TC(kJ/mol)'])
                                                                                                                   
db_nucleotides_cond_reaction['ΔzpE_reaction_d(kJ/mol)'] = (db_nucleotides_cond_reaction['zpE_water(kJ/mol)'] + db_nucleotides_cond_reaction['zpE_nucleotide(kJ/mol)']) - (db_nucleotides_cond_reaction['zpE_backbone(kJ/mol)']+ db_nucleotides_cond_reaction['zpE_RU(kJ/mol)'])
                                                                                                                                                                      
db_nucleotides_cond_reaction['ΔE_reaction(kJ/mol)'] = db_nucleotides_cond_reaction['ΔzpE_reaction_c(kJ/mol)'] + db_nucleotides_cond_reaction['ΔzpE_reaction_d(kJ/mol)']


# For corrE
db_nucleotides_cond_reaction['ΔcorrE_reaction_c(kJ/mol)'] = (db_nucleotides_cond_reaction['corrE_water(kJ/mol)'] + db_nucleotides_cond_reaction['corrE_backbone(kJ/mol)']) - (db_nucleotides_cond_reaction['corrE_IL(kJ/mol)'] + db_nucleotides_cond_reaction['corrE_TC(kJ/mol)'])
                                                          
db_nucleotides_cond_reaction['ΔcorrE_reaction_d(kJ/mol)'] = (db_nucleotides_cond_reaction['corrE_water(kJ/mol)'] + db_nucleotides_cond_reaction['corrE_nucleotide(kJ/mol)']) - (db_nucleotides_cond_reaction['corrE_backbone(kJ/mol)']+ db_nucleotides_cond_reaction['corrE_RU(kJ/mol)'])
                                                                                                              
db_nucleotides_cond_reaction['ΔE(ZPE)_reaction(kJ/mol)'] = db_nucleotides_cond_reaction['ΔcorrE_reaction_c(kJ/mol)'] + db_nucleotides_cond_reaction['ΔcorrE_reaction_d(kJ/mol)']


# For dH
db_nucleotides_cond_reaction['ΔH_reaction_c(kJ/mol)'] = (db_nucleotides_cond_reaction['H_water(kJ/mol)'] + db_nucleotides_cond_reaction['H_backbone(kJ/mol)']) - (db_nucleotides_cond_reaction['H_IL(kJ/mol)'] + db_nucleotides_cond_reaction['H_TC(kJ/mol)'])
                                                          
db_nucleotides_cond_reaction['ΔH_reaction_d(kJ/mol)'] = (db_nucleotides_cond_reaction['H_water(kJ/mol)'] + db_nucleotides_cond_reaction['H_nucleotide(kJ/mol)']) - (db_nucleotides_cond_reaction['H_backbone(kJ/mol)'] + db_nucleotides_cond_reaction['H_RU(kJ/mol)'])
                                                                                                              
db_nucleotides_cond_reaction['ΔH°_reaction(kJ/mol)'] = db_nucleotides_cond_reaction['ΔH_reaction_c(kJ/mol)'] + db_nucleotides_cond_reaction['ΔH_reaction_d(kJ/mol)']


# For TdS
db_nucleotides_cond_reaction['Δ(TdS)_reaction_c(kJ/mol)'] = (db_nucleotides_cond_reaction['TdS_water(kJ/mol)'] + db_nucleotides_cond_reaction['TdS_backbone(kJ/mol)']) - (db_nucleotides_cond_reaction['TdS_IL(kJ/mol)'] + db_nucleotides_cond_reaction['TdS_TC(kJ/mol)'])
                                                            
db_nucleotides_cond_reaction['Δ(TdS)_reaction_d(kJ/mol)'] = (db_nucleotides_cond_reaction['TdS_water(kJ/mol)'] + db_nucleotides_cond_reaction['TdS_nucleotide(kJ/mol)']) - (db_nucleotides_cond_reaction['TdS_backbone(kJ/mol)'] + db_nucleotides_cond_reaction['TdS_RU(kJ/mol)'])
                                                            
db_nucleotides_cond_reaction['Δ(TΔS°)_reaction(kJ/mol)'] = db_nucleotides_cond_reaction['Δ(TdS)_reaction_c(kJ/mol)'] + db_nucleotides_cond_reaction['Δ(TdS)_reaction_d(kJ/mol)']
                                                                                                                    

# For dG
db_nucleotides_cond_reaction['ΔG_reaction_c(kJ/mol)'] = (db_nucleotides_cond_reaction['G_water(kJ/mol)'] + db_nucleotides_cond_reaction['G_backbone(kJ/mol)']) - (db_nucleotides_cond_reaction['G_IL(kJ/mol)'] + db_nucleotides_cond_reaction['G_TC(kJ/mol)'])
                                                          
db_nucleotides_cond_reaction['ΔG_reaction_d(kJ/mol)'] = (db_nucleotides_cond_reaction['G_water(kJ/mol)'] + db_nucleotides_cond_reaction['G_nucleotide(kJ/mol)']) - (db_nucleotides_cond_reaction['G_backbone(kJ/mol)'] + db_nucleotides_cond_reaction['G_RU(kJ/mol)'])
                                                                                                              
db_nucleotides_cond_reaction['ΔG°_reaction(kJ/mol)'] = db_nucleotides_cond_reaction['ΔG_reaction_c(kJ/mol)'] + db_nucleotides_cond_reaction['ΔG_reaction_d(kJ/mol)']  

# Convert pandas dataframe to csv file
db_nucleotides_cond_reaction.to_csv('condensation_reaction_allData_nucleotides_c+d_'+str(IL_name)+'.csv', index=False)                                           

# Read the csv file with all info
db = pd.read_csv ('condensation_reaction_allData_nucleotides_c+d_'+str(IL_name)+'.csv')

# Change name of RUs in the field of this final dataframe
for i, RU in enumerate(RUs_names):
    # Change name of RUs in the field of this final dataframe
    db['RU'] = db['RU'].str.replace(RU,new_RUs_names[i])
    
# Select columns that include only thermo param
for thermo in thermo_params:
    db_to_work = db[["Environment", "TC", "Initial_sugar_ring_conf", "RU", str(thermo)+"_reaction(kJ/mol)"]]
    # Select from previous db columns for specific TC
    for TC in TC_names:
        TC_df = db_to_work[db_to_work['TC'] == TC]
        # Work on furanose sugars
        if TC in ['2deoxy_ribofuranose', 'ribofuranose', 'threose']: # FURANOSES
            # Generate final database TC, RU and env fields
            env_list = np.repeat(np.array(["vacuum", "water"]), [11], axis=0).tolist()
            RUs_list = new_RUs_names + new_RUs_names
            final_db = pd.DataFrame({'Env': env_list, 'TCs': TC, 'RUs': RUs_list})
            # Output entries for specific ring conformation
            for i, ring_conf in enumerate(Ring_Conf_furanose): 
                db_ring_conf = TC_df[TC_df['Initial_sugar_ring_conf'].str.contains(ring_conf)]
                # Add previous field to final db
                final_db[str(new_ring_conf_fura[i])] = db_ring_conf[str(thermo)+'_reaction(kJ/mol)'].to_numpy()
            # Obtain freq hist
            gen_bar_graph (final_db, thermo, TC, IL_name)
            final_db.to_csv('summary_'+str(thermo)+'_'+str(TC)+'.csv', index=False)
        elif TC in ["2deoxy_ribopyranose", "ribopyranose"]: ## PYRANOSES   
            # Generate final database TC, RU and env fields
            env_list = np.repeat(np.array(["vacuum", "water"]), [11], axis=0).tolist()
            RUs_list = new_RUs_names + new_RUs_names
            final_db = pd.DataFrame({'Env': env_list, 'TCs': TC, 'RUs': RUs_list})            
            # Output entries for specific ring conformation
            for i, ring_conf in enumerate(Ring_Conf_pyranose): 
                db_ring_conf = TC_df[TC_df['Initial_sugar_ring_conf'].str.contains(ring_conf)]
                # Add previous field to final db
                final_db[str(new_ring_conf_pyra[i])] = db_ring_conf[str(thermo)+'_reaction(kJ/mol)'].to_numpy()
            # Obtain freq hist
            gen_bar_graph (final_db, thermo, TC, IL_name)
            final_db.to_csv('summary_'+str(thermo)+'_'+str(TC)+'.csv', index=False)        
        elif TC in ["glycerol", "glyceric_acid"]: ## NON SUGAR            
            # Generate final database TC, RU and env fields
            env_list = np.repeat(np.array(["vacuum", "water"]), [11], axis=0).tolist()
            RUs_list = new_RUs_names + new_RUs_names
            final_db = pd.DataFrame({'Env': env_list, 'TCs': TC, 'RUs': RUs_list})       
            db_ring_conf = TC_df[TC_df['Initial_sugar_ring_conf'] == 'none']           
            # Add previous field to final db
            final_db[str(thermo)+'_react'] = db_ring_conf[str(thermo)+'_reaction(kJ/mol)'].to_numpy()            
            # Obtain freq hist     
            gen_bar_graph (final_db, thermo, TC, IL_name)
            final_db.to_csv('summary_'+str(thermo)+'_'+str(TC)+'.csv', index=False)
            
# Move files to their new folders
for TC in TC_names:
    os.mkdir(TC)
    for thermo in thermo_params:
        shutil.move('Condensation_react_'+str(TC)+'_'+str(thermo)+'_'+str(IL_name)+'.pdf', str(TC)+'/Condensation_react_'+str(TC)+'_'+str(thermo)+'_'+str(IL_name)+'.pdf')
        shutil.move('summary_'+str(thermo)+'_'+str(TC)+'.csv', str(TC)+'/summary_'+str(thermo)+'_'+str(TC)+'.csv')
            
END_SCRIPT

# Run script.py with python3
python3 create_csv_files_hist_cond_reaction.py 

cd ../ # LEAVE IL FOLDER

done














