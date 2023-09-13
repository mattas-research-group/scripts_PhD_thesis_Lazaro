#!/bin/bash

#       f_pseudorot_equil_nucleosides.sh

#	09-jun-2023	Fixed function calc_deltas in gen_pseudo_rot_equil_csv_hist.py
# 	23-Jan-2023 	Created

# Description:
# This script obtains the delta thermo params for the 
# pseudorotation equilibrium of sugars rings.

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a RU_names=("adenine" "guanine" "cytosine" "thymine" "uracil" "TARC_cbond" "TARC_nbond" "BA_cbond" "BA_nbond" "CA" "melamine")


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

# Load function input variables
local environment=$1
local TC=$2
local ring_conf=$3
local RU=$4

nucleoside_name="$TC"_"$ring_conf"_"$RU"

############################################# OBTAIN ZERO POINT ENERGY ################################

# Assign variable to zero point energy in kJ/mol txt file
zero_point_energy=`cat "$nucleoside_name"_zeropoint_energy_kJ_mol.txt`

############################################# OBTAIN CORRECTED ENERGY ################################

# Assign variable to zero point energy in kJ/mol txt file
corrected_energy=`cat "$nucleoside_name"_corrected_energy_kJ_mol.txt`

############################################# OBTAIN ENTHALPY ################################

# Assign variable to zero point energy in kJ/mol txt file
enthalpy=`cat "$nucleoside_name"_entalphy_kJ_mol.txt`

############################################ OBTAIN THE FREE ENERGY ##################################

# Assign variable to zero point energy in kJ/mol txt file
free_energy=`cat "$nucleoside_name"_G_kJ_mol.txt`

############################################ OBTAIN PARAMETERS FOR WATER ##############################

# OUTPUT ENERGIES AND OTHER INFO IN CORRESPONDING SUMMARY FILES IN ANALYSIS FOLDER
if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then
echo $zero_point_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$TC/"$ring_conf_fura"_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$TC/"$ring_conf_fura"_corrected_energies.txt
echo $enthalpy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$TC/"$ring_conf_fura"_enthalpies.txt
echo $free_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$TC/"$ring_conf_fura"_free_energies.txt
else
echo $zero_point_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$TC/"$ring_conf_pyra"_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$TC/"$ring_conf_pyra"_corrected_energies.txt
echo $enthalpy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$TC/"$ring_conf_pyra"_enthalpies.txt
echo $free_energy >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$TC/"$ring_conf_pyra"_free_energies.txt
fi

#################################################################################### END OF FUNCTION ###########################################################################################

}

########################################################
##                                                    ##
## FUNCTION TO ORDER AND CREATE FINAL CSV FILES WITH  ##
## ALL THERMODYNAMIC FOR RING CONFIG EQUILIBRIUM.     ##
##                                                    ##
########################################################

function create_csv_RingConf_equil {

local TC=$1

cat >create_csv_files_for_RingConf_equil.py <<'END_SCRIPT'
# Import environmental variables from python libraries
import sys, os
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from numpy import loadtxt
import seaborn as sb
import scipy as scipy
from scipy import stats



END_SCRIPT

# Run script.py with python3
#python3 create_csv_files_for_RingConf_equil.py



}



##########################################################################################################################
############################################## WORK ON EACH FOLDER #######################################################
##########################################################################################################################

# Create folder for analysis and post-processing
mkdir $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium


# Open corresponding folder
cd nucleosides/final_gaussian_opt

for environment in "${environment[@]}"; do

# Create folder in post_proc_folder
#mkdir $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$environment


cd $environment

for TC in "${TC_names[@]}"; do

cd $TC

if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

# Create corresponding folder for furanoses
mkdir $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$TC

for ring_conf_fura in "${Ring_Conf_furanose[@]}"; do

cd $ring_conf_fura

for RU in "${RU_names[@]}"; do

cd $RU

# Output message of status of script
echo "working on $environment/$TC/$ring_conf_fura/$RU"

# Run function to extract thermo parameters from log files for nucleoside, RU and TC 
extract_thermo_params_pseu_rot $environment $TC $ring_conf_fura $RU 


cd ../
done

cd ../
done

elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then

# Create corresponding folder for furanoses
mkdir $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$TC


for ring_conf_pyra in "${Ring_Conf_pyranose[@]}"; do

cd $ring_conf_pyra

for RU in "${RU_names[@]}"; do

cd $RU

# Output message of status of script
echo "working on $environment/$TC/$ring_conf_pyra/$RU"

# Run function to extract thermo parameters from log files for nucleoside, RU and TC 
extract_thermo_params_pseu_rot $environment $TC $ring_conf_pyra $RU 



cd ../
done

cd ../
done

fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT

cd ../ # LEAVE FINAL_GAUSSIAN_OPT FOLDER AND NUCLEOSIDES FOLDER

#################################################################################################################################################
############################################## Generate final csv files with equil params #######################################################
#################################################################################################################################################

cd post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium

#################################################################################################################################################
############################################## Generate final csv files with equil params #######################################################
#################################################################################################################################################
for TC in "${TC_names[@]}"; do

cd $TC

export TC

# Create and run python script to gen summary csv files and freq hist 
cat >gen_pseudo_rot_equil_csv_hist.py <<'END_SCRIPT'
import sys, os
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from numpy import loadtxt
import seaborn as sb
import scipy as scipy
from scipy import stats
from matplotlib.gridspec import SubplotSpec

################################### Define necesary variables (stacked lists) ###################################
environment = ['vacuum','water']
TC_names = ["2deoxy_ribopyranose", "ribopyranose"]
RUs_names = ["adenine", "guanine", "cytosine", "thymine", "uracil", "TARC_cbond", "TARC_nbond", "BA_cbond", "BA_nbond", "CA", "melamine"]
new_RUs_names = ["A", "G", "C", "T", "U", "cTARC", "nTARC", "cBA", "nBA", "CA", "MM"]
type_RU = ["canonical", "canonical", "canonical"]
Ring_Conf_furanose = ["2endo_alpha", "2endo_beta", "3endo_alpha", "3endo_beta"]
Ring_Conf_pyranose = ["D1C4_alpha", "D1C4_beta", "D4C1_alpha", "D4C1_beta"]

# Import environmental variables from bash script for angles in degrees
TC = os.environ.get("TC")

################################## Generate necesary functions #############################################

########################################################################
##                                                                    ##
## FUNCTIONS TO GENERATE DELTA THERMO PARAMS.                         ##
##                                                                    ##
## Description:                                                       ##
## This function generates for each thermo                            ##
## param and each TC delta values for each                            ##
## pseudorot equil and anomer-exchange react.                         ##
##                                                                    ##
## Input                                                              ##
## zpE		: single point energy (0K)                             ##
## corrE	: corrected energy (298K)                              ##
## dH		: enthalpy                                             ##
## dG		: Gibbs energy                                         ##
## TdS		: T * dEntropy                                         ##
##                                                                    ##
## for the following rings:                                           ##
## fura:  ["2endo_alpha", "2endo_beta", "3endo_alpha", "3endo_beta"]  ##
## pyra:  ["D1C4_alpha", "D1C4_beta", "D4C1_alpha", "D4C1_beta"]      ##
##                                                                    ##
########################################################################

def calc_deltas(zpE, corrE, dH, dG, TdS, TC):
    thermo_v = [zpE, corrE, dH, dG, TdS]

    zpE_deltas = [] 
    corrE_deltas = [] 
    dH_deltas = [] 
    dG_deltas = [] 
    TdS_deltas = []
    
    deltas = [zpE_deltas, corrE_deltas, dH_deltas, dG_deltas, TdS_deltas] # [[], [], [], [], []]

    for i, thermo_param in enumerate(thermo_v): 
        if TC == "2deoxy_ribofuranose" or TC == "ribofuranose" or TC == "threose":
            # ["2endo_alpha", "2endo_beta", "3endo_alpha", "3endo_beta"]  
            d1 = thermo_param[0] - thermo_param[2] 
            d2 = thermo_param[1] - thermo_param[0] 
            d3 = thermo_param[1] - thermo_param[3] 
            d4 = thermo_param[3] - thermo_param[2] 
            d5 = thermo_param[1] - thermo_param[2] 
            d6 = thermo_param[0] - thermo_param[3]
        else:  
	    # ["D1C4_alpha", "D1C4_beta", "D4C1_alpha", "D4C1_beta"]
            d1 = thermo_param[2] - thermo_param[0] 
            d2 = thermo_param[3] - thermo_param[2] 
            d3 = thermo_param[3] - thermo_param[1] 
            d4 = thermo_param[1] - thermo_param[0] 
            d5 = thermo_param[3] - thermo_param[0] 
            d6 = thermo_param[2] - thermo_param[1]
		    
        for d in [d1, d2, d3, d4, d5, d6]:
            deltas[i].append(d)

    return deltas             


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
## param for pseudorotation equil and anomer    ##
## exchange reactions (6 reactions in total).   ##                   
##                                              ##
##################################################
  
def gen_freq_hist (thermo):

    # Define number of rows and cols
    rows = 2
    cols = 2

    # Create the fig and axn
    fig, axn = plt.subplots(rows, cols, sharex=False, sharey=False, figsize=(16,10)) 
    
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
    
    for i, (ax, env, type_ru) in enumerate([(axn.flat[0], "vacuum", "canonical"), 
                                            (axn.flat[1], "vacuum", "non_canonical"), 
                                            (axn.flat[2], "water", "canonical"),
                                            (axn.flat[3], "water", "non_canonical")]):       
        # Obtain rows from database that belong to vacuum or water env
        env_df = dataframe[dataframe['Env:'].isin([env])]    
        RUs_canonical = ["A", "G", "C", "T", "U"]
        RUs_non_canonical = ["cTARC", "nTARC", "cBA", "nBA", "CA", "MM"]
        # From this new df obtain rows for either canonical or non canonical bases
        if type_ru == "canonical":
            ru_df = env_df[env_df['RU_name'].isin(RUs_canonical)]
        else:
            ru_df = env_df[env_df['RU_name'].isin(RUs_non_canonical)]        
        # Set indexes of database
        ru_df = ru_df.set_index("RU_name")
        #print(ru_df)    
        # Generate from this df the columns for dG values + RU + environment
        dG_list = []
        for i in range(1, 7, 1):
            G = str(thermo)+"#"+str(i)
            dG_list.append(G)
        dG_df = ru_df.filter(dG_list, axis=1)   
        # Generate bar graph now
        if env == "vacuum" and type_ru == "canonical":
            ax1 = dG_df.plot(kind='bar', ax=ax, colormap='Greys', rot=360, width=0.8, edgecolor = "black", 
                             xlabel='', ylabel='ΔG (kJ/mol)', legend=False, fontsize=20)
            ax1.set_xticks([]) 
            ax1.set_ylabel(str(thermo)+"$^\\circ$ (kJ/mol)", fontsize=20, fontdict=dict(weight='bold'))
            ax1.set_title('Canonical RUs', weight='bold', fontsize=22)
            for container in ax1.containers:
                ax1.bar_label(container, fmt='%.1f', fontsize=14)
        elif env == "vacuum" and type_ru == "non_canonical":
            ax2 = dG_df.plot(kind='bar', ax=ax, colormap='Greys', width=0.8, rot=360, edgecolor = "black", 
                             xlabel='', legend=False, title='Vacuum (non_canonical RUs)', fontsize=20)
            ax2.set_xticks([])
            ax2.set_title('Non-canonical RUs', weight='bold', fontsize=22)
            for container in ax2.containers:
                ax2.bar_label(container, fmt='%.1f', fontsize=14)
        elif env == "water" and type_ru == "canonical":
            ax3 = dG_df.plot(kind='bar', ax=ax, colormap='Greys', width=0.8, rot=360, edgecolor = "black", 
                             xlabel='', ylabel='ΔG (kJ/mol)', legend=False, fontsize=20)
            ax3.set_ylabel(str(thermo) + "$^\\circ$ (kJ/mol)", fontsize=20, fontdict=dict(weight='bold'))
            #ax3.set_title('Water canonical RUs', weight='bold', fontsize=18)
            for container in ax3.containers:
                ax3.bar_label(container, fmt='%.1f', fontsize=14)
        elif env == "water" and type_ru == "non_canonical":
            ax4 = dG_df.plot(kind='bar', ax=ax, colormap='Greys', width=0.8, rot=360, edgecolor = "black", 
                             xlabel='', legend=True, fontsize=20)
            #ax4.set_title('Water non_canonical RUs', weight='bold', fontsize=18)
            ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=18)
            for container in ax4.containers:
                ax4.bar_label(container, fmt='%.1f', fontsize=16)            
    # Set parameters for axis labels and titles       
    for axis in [ax1, ax2, ax3, ax4]:
        axis.axhline(y = 0.0, color = 'black', linestyle = '--')
        labels = axis.get_xticklabels() 
        #for label in labels:
            #label.set_fontweight('bold')    
    plt.tight_layout()
    plt.savefig('pseudo_rot_'+str(thermo)+'.pdf')
    
zpE = []
corrE = []
dH = []
dG = []
TdS = []

if TC == "2deoxy_ribofuranose" or TC == "ribofuranose" or TC == "threose":
    for fura_Rconf in Ring_Conf_furanose:        
        zpE.append(loadtxt(str(fura_Rconf) + '_zero_point_energies.txt', dtype=float)) # Append to general list
        corrE.append(loadtxt(str(fura_Rconf) + '_corrected_energies.txt', dtype=float)) # Append to general list               
        dH.append(loadtxt(str(fura_Rconf) + '_enthalpies.txt', dtype=float)) # Append to general list                             
        dG.append(loadtxt(str(fura_Rconf) + '_free_energies.txt', dtype=float)) # Append to general list               
        TdS.append((loadtxt(str(fura_Rconf) + '_enthalpies.txt', dtype=float)) - (loadtxt(str(fura_Rconf) + '_free_energies.txt', dtype=float)))
elif TC == "2deoxy_ribopyranose" or TC == "ribopyranose":
    for pyra_Rconf in Ring_Conf_pyranose:
        zpE.append(loadtxt(str(pyra_Rconf) + '_zero_point_energies.txt', dtype=float)) # Append to general list
        corrE.append(loadtxt(str(pyra_Rconf) + '_corrected_energies.txt', dtype=float)) # Append to general list               
        dH.append(loadtxt(str(pyra_Rconf) + '_enthalpies.txt', dtype=float)) # Append to general list                           
        dG.append(loadtxt(str(pyra_Rconf) + '_free_energies.txt', dtype=float)) # Append to general list               
        TdS.append((loadtxt(str(pyra_Rconf) + '_enthalpies.txt', dtype=float)) - (loadtxt(str(pyra_Rconf) + '_free_energies.txt', dtype=float)))
        
[zpE_deltas, corrE_deltas, dH_deltas, dG_deltas, TdS_deltas] = calc_deltas (zpE, corrE, dH, dG, TdS, TC)  
        
# Generate csv file
TC_list = [TC] * 22 #len(RUs_names)
env_list = sorted(environment * 11) #len(RUs_names)
RUs_list = new_RUs_names * 2
        
if TC == "2deoxy_ribofuranose" or TC == "ribofuranose" or TC == "threose":
    dataframe = pd.DataFrame({'Env:':env_list, 'TC_name:': TC_list,'RU_name': RUs_list, 
                              'zpE_2endo_alpha': zpE[0], 'zpE_2endo_beta': zpE[1], 
                              'zpE_3endo_alpha': zpE[2], 'zpE_3endo_beta': zpE[3], 
	                      'corrE_2endo_alpha': corrE[0], 'corrE_2endo_beta': corrE[1], 
	                      'corrE_3endo_alpha': corrE[2], 'corrE_3endo_beta': corrE[3], 
	                      'dH_2endo_alpha': dH[0], 'dH_2endo_beta': dH[1], 'dH_3endo_alpha': dH[2], 'dH_3endo_beta': dH[3],
	                      'dG_2endo_alpha': dG[0], 'dG_2endo_beta': dG[1], 'dG_3endo_alpha': dG[2], 'dG_3endo_beta': dG[3], 
	                      'TdS_2endo_alpha': TdS[0], 'TdS_2endo_beta': TdS[1], 'TdS_3endo_alpha': TdS[2], 'TdS_3endo_beta': TdS[3],
	                      'zpΔE#1': zpE_deltas[0], 'zpΔE#2': zpE_deltas[1], 'zpΔE#3': zpE_deltas[2], 
	                      'zpΔE#4': zpE_deltas[3], 'zpΔE#5': zpE_deltas[4], 'zpΔE#6': zpE_deltas[5],
	                      'corrΔE#1': corrE_deltas[0], 'corrΔE#2': corrE_deltas[1], 'corrΔE#3': corrE_deltas[2], 
	                      'corrΔE#4': corrE_deltas[3], 'corrΔE#5': corrE_deltas[4], 'corrΔE#6': corrE_deltas[5],
	                      'ΔH#1': dH_deltas[0], 'ΔH#2': dH_deltas[1], 'ΔH#3': dH_deltas[2], 'ΔH#4': dH_deltas[3], 'ΔH#5': dH_deltas[4], 'ΔH#6': dH_deltas[5],
	                      'ΔG#1': dG_deltas[0], 'ΔG#2': dG_deltas[1], 'ΔG#3': dG_deltas[2], 'ΔG#4': dG_deltas[3], 'ΔG#5': dG_deltas[4], 'ΔG#6': dG_deltas[5],
	                      'Δ(TΔS)#1': TdS_deltas[0], 'Δ(TΔS)#2': TdS_deltas[1], 'Δ(TΔS)#3': TdS_deltas[2], 'Δ(TΔS)#4': TdS_deltas[3], 'Δ(TΔS)#5': TdS_deltas[4], 'Δ(TΔS)#6': TdS_deltas[5]})
	             
    dataframe.to_csv("thermo_fura_"+str(TC)+"_pseudo_eq_allData.csv", index=False) 

elif TC == "2deoxy_ribopyranose" or TC == "ribopyranose":  
    dataframe = pd.DataFrame({'Env:':env_list, 'TC_name:': TC_list,'RU_name': RUs_list, 
                              'zpE_D1C4_alpha': zpE[0], 'zpE_D1C4_beta': zpE[1], 'zpE_D4C1_alpha': zpE[2], 'zpE_D4C1_beta': zpE[3], 
	                       'corrE_D1C4_alpha': corrE[0], 'corrE_D1C4_beta': corrE[1], 'corrE_D4C1_alpha': corrE[2], 'corrE_D4C1_beta': corrE[3], 
	                       'dH_D1C4_alpha': dH[0], 'dH_D1C4_beta': dH[1], 'dH_D4C1_alpha': dH[2], 'dH_D4C1_beta': dH[3],
	                       'dG_D1C4_alpha': dG[0], 'dG_D1C4_beta': dG[1], 'dG_D4C1_alpha': dG[2], 'dG_D4C1_beta': dG[3], 
	                       'TdS_D1C4_alpha': TdS[0], 'TdS_D1C4_beta': TdS[1], 'TdS_D4C1_alpha': TdS[2], 'TdS_D4C1_beta': TdS[3],
	                       'zpΔE#1': zpE_deltas[0], 'zpΔE#2': zpE_deltas[1], 'zpΔE#3': zpE_deltas[2], 'zpΔE#4': zpE_deltas[3], 'zpΔE#5': zpE_deltas[4], 'zpΔE#6': zpE_deltas[5],
	                       'corrΔE#1': corrE_deltas[0], 'corrΔE#2': corrE_deltas[1], 'corrΔE#3': corrE_deltas[2], 'corrΔE#4': corrE_deltas[3], 'corrΔE#5': corrE_deltas[4], 
                               'corrΔE#6': corrE_deltas[5],
	                       'ΔH#1': dH_deltas[0], 'ΔH#2': dH_deltas[1], 'ΔH#3': dH_deltas[2], 'ΔH#4': dH_deltas[3], 'ΔH#5': dH_deltas[4], 'ΔH#6': dH_deltas[5],
	                       'ΔG#1': dG_deltas[0], 'ΔG#2': dG_deltas[1], 'ΔG#3': dG_deltas[2], 'ΔG#4': dG_deltas[3], 'ΔG#5': dG_deltas[4], 'ΔG#6': dG_deltas[5],
	                       'Δ(TΔS)#1': TdS_deltas[0], 'Δ(TΔS)#2': TdS_deltas[1], 'Δ(TΔS)#3': TdS_deltas[2], 'Δ(TΔS)#4': TdS_deltas[3], 'Δ(TΔS)#5': TdS_deltas[4], 'Δ(TΔS)#6': TdS_deltas[5]})
                    
    dataframe.to_csv("thermo_pyra_"+str(TC)+"_pseudo_eq_allData.csv", index=False) 
        
# Generate the freq hist plots
thermo_params = ['zpΔE', 'corrΔE', 'ΔH', 'Δ(TΔS)', 'ΔG']
for thermo in thermo_params:
    gen_freq_hist (thermo)     



END_SCRIPT

# Run script.py with python3
python3 gen_pseudo_rot_equil_csv_hist.py

cd ../

done











