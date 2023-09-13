#!/bin/bash


#	f_pseudorot_equil_nucleotides_c+d.sh

#	18-jun-2023	Fixed function calc_deltas in gen_pseudo_rot_equil_csv_hist.py
# 	04-mar-2023 	Created

# Description:
# This script obtains the delta thermo params for the 
# pseudorotation equilibrium of sugars rings.

# To run script type:
# ./f_pseudorot_equil_nucleotides_c+d.sh

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a RU_names=("adenine" "guanine" "cytosine" "thymine" "uracil" "TARC_cbond" "TARC_nbond" "BA_cbond" "BA_nbond" "CA" "melamine")
declare -a IL_names=("phosphate" "arsenate")

##################################################
##                                              ##
## FUNCTIONS TO OBTAIN DIFFERENT THERMODYNAMIC  ##
## PARAMETERS FOR EACH RING_CONF OF TC(SUGARS)  ##
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
local IL=$5

# Declare nucleotide name
nucleotide_name="$TC"_"$ring_conf"_"$RU"_"$IL"

############################################# OBTAIN ZERO POINT ENERGY ################################

# Assign variable to zero point energy in kcal/mol txt file
zero_point_energy=`cat "$nucleotide_name"_zeropoint_energy_kj_mol.txt`

############################################# OBTAIN CORRECTED ENERGY ################################

# Assign variable to zero point energy in kcal/mol txt file
corrected_energy=`cat "$nucleotide_name"_corrected_energy_kj_mol.txt`

############################################# OBTAIN ENTHALPY ################################

# Assign variable to zero point energy in kcal/mol txt file
enthalpy=`cat "$nucleotide_name"_entalphy_kj_mol.txt`

############################################ OBTAIN THE FREE ENERGY ##################################

# Assign variable to zero point energy in kcal/mol txt file
free_energy=`cat "$nucleotide_name"_G_kj_mol.txt`

# OUTPUT ENERGIES AND OTHER INFO IN CORRESPONDING SUMMARY FILES IN ANALYSIS FOLDER
echo $zero_point_energy >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$IL/$TC/"$ring_conf"_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$IL/$TC/"$ring_conf"_corrected_energies.txt
echo $enthalpy >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$IL/$TC/"$ring_conf"_enthalpies.txt
echo $free_energy >> $MD_RUN/nucleotides_c+d/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$IL/$TC/"$ring_conf"_free_energies.txt


#################################################################################### END OF FUNCTION ###########################################################################################

}







##########################################################################################################################
############################################## WORK ON EACH FOLDER #######################################################
##########################################################################################################################

# Open corresponding folder
mkdir nucleotides_c+d/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium

# Create IL folder inside the post-processing folder
for IL in "${IL_names[@]}"; do
mkdir nucleotides_c+d/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$IL

# Create TC folder inside the post-processing folder
for TC in "${TC_names[@]}"; do
mkdir nucleotides_c+d/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/$IL/$TC
done

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

# Run function to extract thermo parameters from log files for nucleoside, RU and TC 
extract_thermo_params_pseu_rot $environment $TC $ring_conf_fura $RU $IL





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

# Run function to extract thermo parameters from log files for nucleoside, RU and TC 
extract_thermo_params_pseu_rot $environment $TC $ring_conf_pyra $RU $IL







cd ../
done

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

#################################################################################################################################################
############################################## Generate final csv files with equil params #######################################################
#################################################################################################################################################

# Open folder
cd post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium

for IL in "${IL_names[@]}"; do

cd $IL

for TC in "${TC_names[@]}"; do

cd $TC

export IL
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

environment = ['vacuum','water']
TC_names = ["2deoxy_ribopyranose", "ribopyranose"]
RUs_names = ["adenine", "guanine", "cytosine", "thymine", "uracil", "TARC_cbond", "TARC_nbond", "BA_cbond", "BA_nbond", "CA", "melamine"]
new_RUs_names = ["A", "G", "C", "T", "U", "TAP-$C^5$", "TAP-N", "BA-$C^5$", "BA-N", "CA", "MM"]
type_RU = ["canonical", "canonical", "canonical"]
Ring_Conf_furanose = ["2endo_alpha", "2endo_beta", "3endo_alpha", "3endo_beta"]
Ring_Conf_pyranose = ["D1C4_alpha", "D1C4_beta", "D4C1_alpha", "D4C1_beta"]

# Import environmental variables from bash script for angles in degrees
TC = os.environ.get("TC")
IL = os.environ.get("IL")

# Create necesary functions
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

def calc_deltas(zpE, corrE, dH, dG, TdS):
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

def gen_freq_hist (dataframe, thermo, TC, IL):
    fig, axn = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(22,12)) 
    for i, (ax, env, type_ru) in enumerate([(axn.flat[0], "vacuum", "canonical"), 
                                            (axn.flat[1], "vacuum", "non_canonical"), 
                                            (axn.flat[2], "water", "canonical"),
                                            (axn.flat[3], "water", "non_canonical")]):      
                                             
        # Obtain rows from database that belong to vacuum or water env
        env_df = dataframe[dataframe['Env:'].isin([env])]    
        RUs_canonical = ["A", "G", "C", "T", "U"]
        RUs_non_canonical = ["TAP-$C^5$", "TAP-N", "BA-$C^5$", "BA-N", "CA", "MM"]
        
        # From this new df obtain rows for either canonical or non canonical bases
        if type_ru == "canonical":
            ru_df = env_df[env_df['RU_name'].isin(RUs_canonical)]
        else:
            ru_df = env_df[env_df['RU_name'].isin(RUs_non_canonical)]        
        # Set indexes of database
        ru_df = ru_df.set_index("RU_name")
        
        # Generate from this df the columns for dG values + RU + environment
        dG_list = []
        for i in range(1, 7, 1):
            if thermo == "ΔE(ZPE)":
                G = "ΔE$_"+str(i)+"(ZPE)$"
                dG_list.append(G)
            else:
                G = str(thermo)+"$_"+str(i)+"$"
                dG_list.append(G) 
        dG_df = ru_df.filter(dG_list, axis=1) 
          
        # Generate bar graph now
        if env == "vacuum" and type_ru == "canonical":
            ax1 = dG_df.plot(kind='bar', ax=ax, colormap='Greys', rot=360, width=0.8, edgecolor = "black", 
                             xlabel='', ylabel='ΔG (kJ/mol)', legend=False, fontsize=28)
            ax1.set_xticks([]) 
            ax1.set_ylabel(str(thermo)+"(kJ/mol)", fontsize=28, fontdict=dict(weight='bold'))
            ax1.set_title('Canonical RUs', weight='bold', fontsize=32)
            for container in ax1.containers:
                ax1.bar_label(container, fmt='%.1f', fontsize=18)
            # Generate vertical lines
            [ax1.axvline(x, color = 'grey', linestyle='--') for x in [0.5, 1.5, 2.5, 3.5]]
        elif env == "vacuum" and type_ru == "non_canonical":
            ax2 = dG_df.plot(kind='bar', ax=ax, colormap='Greys', width=0.8, rot=360, edgecolor = "black", 
                             xlabel='', legend=False, title='Vacuum (non_canonical RUs)', fontsize=28)
            ax2.set_xticks([])
            ax2.set_title('Non-canonical RUs', weight='bold', fontsize=32)
            for container in ax2.containers:
                ax2.bar_label(container, fmt='%.1f', fontsize=18)
            # Generate vertical lines
            [ax2.axvline(x, color = 'grey', linestyle='--') for x in [0.5, 1.5, 2.5, 3.5, 4.5]] 
        elif env == "water" and type_ru == "canonical":
            ax3 = dG_df.plot(kind='bar', ax=ax, colormap='Greys', width=0.8, rot=360, edgecolor = "black", 
                             xlabel='', ylabel='ΔG (kJ/mol)', legend=False, fontsize=28)
            ax3.set_ylabel(str(thermo) + "(kJ/mol)", fontsize=28, fontdict=dict(weight='bold'))
            #ax3.set_title('Water canonical RUs', weight='bold', fontsize=18)
            for container in ax3.containers:
                ax3.bar_label(container, fmt='%.1f', fontsize=18)
            # Generate vertical lines
            [ax3.axvline(x, color = 'grey', linestyle='--') for x in [0.5, 1.5, 2.5, 3.5]] 
        elif env == "water" and type_ru == "non_canonical":
            ax4 = dG_df.plot(kind='bar', ax=ax, colormap='Greys', width=0.8, rot=360, edgecolor = "black", 
                             xlabel='', legend=True, fontsize=28)
            #ax4.set_title('Water non_canonical RUs', weight='bold', fontsize=18)
            ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=24)
            for container in ax4.containers:
                ax4.bar_label(container, fmt='%.1f', fontsize=18)
            # Generate vertical lines
            [ax4.axvline(x, color = 'grey', linestyle='--') for x in [0.5, 1.5, 2.5, 3.5, 4.5]]                
    # Set parameters for axis labels and titles       
    for axis in [ax1, ax2, ax3, ax4]:
        axis.axhline(y = 0.0, color = 'black', linestyle = '--')
        labels = axis.get_xticklabels() 
        #for label in labels:
            #label.set_fontweight('bold')    
    plt.tight_layout()
    plt.savefig('pseudo_rot_'+str(thermo)+'_'+str(TC)+'_'+str(IL)+'_c+d.pdf')
    
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
        
[zpE_deltas, corrE_deltas, dH_deltas, dG_deltas, TdS_deltas] = calc_deltas (zpE, corrE, dH, dG, TdS)  
        
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
	                      'ΔE$_1$': zpE_deltas[0], 'ΔE$_2$': zpE_deltas[1], 'ΔE$_3$': zpE_deltas[2], 'ΔE$_4$': zpE_deltas[3], 'ΔE$_5$': zpE_deltas[4], 'ΔE$_6$': zpE_deltas[5],
	                      'ΔE$_1(ZPE)$': corrE_deltas[0], 'ΔE$_2(ZPE)$': corrE_deltas[1], 'ΔE$_3(ZPE)$': corrE_deltas[2], 'ΔE$_4(ZPE)$': corrE_deltas[3], 'ΔE$_5(ZPE)$': corrE_deltas[4], 'ΔE$_6(ZPE)$': corrE_deltas[5],
	                      'ΔH°$_1$': dH_deltas[0], 'ΔH°$_2$': dH_deltas[1], 'ΔH°$_3$': dH_deltas[2], 'ΔH°$_4$': dH_deltas[3], 'ΔH°$_5$': dH_deltas[4], 'ΔH°$_6$': dH_deltas[5],
	                      'ΔG°$_1$': dG_deltas[0], 'ΔG°$_2$': dG_deltas[1], 'ΔG°$_3$': dG_deltas[2], 'ΔG°$_4$': dG_deltas[3], 'ΔG°$_5$': dG_deltas[4], 'ΔG°$_6$': dG_deltas[5],
	                      'Δ(TΔS°)$_1$': TdS_deltas[0], 'Δ(TΔS°)$_2$': TdS_deltas[1], 'Δ(TΔS°)$_3$': TdS_deltas[2], 'Δ(TΔS°)$_4$': TdS_deltas[3], 'Δ(TΔS°)$_5$': TdS_deltas[4], 'Δ(TΔS°)$_6$': TdS_deltas[5]})
	             
    dataframe.to_csv("thermo_fura_"+str(TC)+"_pseudo_eq_allData.csv", index=False) 

elif TC == "2deoxy_ribopyranose" or TC == "ribopyranose":  
    dataframe = pd.DataFrame({'Env:':env_list, 'TC_name:': TC_list,'RU_name': RUs_list, 
                              'zpE_D1C4_alpha': zpE[0], 'zpE_D1C4_beta': zpE[1], 'zpE_D4C1_alpha': zpE[2], 'zpE_D4C1_beta': zpE[3], 
	                       'corrE_D1C4_alpha': corrE[0], 'corrE_D1C4_beta': corrE[1], 'corrE_D4C1_alpha': corrE[2], 'corrE_D4C1_beta': corrE[3], 
	                       'dH_D1C4_alpha': dH[0], 'dH_D1C4_beta': dH[1], 'dH_D4C1_alpha': dH[2], 'dH_D4C1_beta': dH[3],
	                       'dG_D1C4_alpha': dG[0], 'dG_D1C4_beta': dG[1], 'dG_D4C1_alpha': dG[2], 'dG_D4C1_beta': dG[3], 
	                       'TdS_D1C4_alpha': TdS[0], 'TdS_D1C4_beta': TdS[1], 'TdS_D4C1_alpha': TdS[2], 'TdS_D4C1_beta': TdS[3],
	                       'ΔE$_1$': zpE_deltas[0], 'ΔE$_2$': zpE_deltas[1], 'ΔE$_3$': zpE_deltas[2], 'ΔE$_4$': zpE_deltas[3], 'ΔE$_5$': zpE_deltas[4], 'ΔE$_6$': zpE_deltas[5],
	                       'ΔE$_1(ZPE)$': corrE_deltas[0], 'ΔE$_2(ZPE)$': corrE_deltas[1], 'ΔE$_3(ZPE)$': corrE_deltas[2], 'ΔE$_4(ZPE)$': corrE_deltas[3], 'ΔE$_5(ZPE)$': corrE_deltas[4], 
                               'ΔE$_6(ZPE)$': corrE_deltas[5],
	                       'ΔH°$_1$': dH_deltas[0], 'ΔH°$_2$': dH_deltas[1], 'ΔH°$_3$': dH_deltas[2], 'ΔH°$_4$': dH_deltas[3], 'ΔH°$_5$': dH_deltas[4], 'ΔH°$_6$': dH_deltas[5],
	                       'ΔG°$_1$': dG_deltas[0], 'ΔG°$_2$': dG_deltas[1], 'ΔG°$_3$': dG_deltas[2], 'ΔG°$_4$': dG_deltas[3], 'ΔG°$_5$': dG_deltas[4], 'ΔG°$_6$': dG_deltas[5],
	                       'Δ(TΔS°)$_1$': TdS_deltas[0], 'Δ(TΔS°)$_2$': TdS_deltas[1], 'Δ(TΔS°)$_3$': TdS_deltas[2], 'Δ(TΔS°)$_4$': TdS_deltas[3], 'Δ(TΔS°)$_5$': TdS_deltas[4], 'Δ(TΔS°)$_6$': TdS_deltas[5]})
                    
    dataframe.to_csv("thermo_pyra_"+str(TC)+"_pseudo_eq_allData.csv", index=False)
        
# Generate the freq hist plots
thermo_params = ['ΔE', 'ΔE(ZPE)', 'ΔH°', 'ΔG°', 'Δ(TΔS°)']
for thermo in thermo_params:
    gen_freq_hist (dataframe, thermo, TC, IL)   
    
################# Merge dE, dE(ZPE) and dG values in one cell for tables in paper  
if TC == "2deoxy_ribofuranose" or TC == "ribofuranose" or TC == "threose":
    df = pd.read_csv('thermo_fura_'+str(TC)+'_pseudo_eq_allData.csv')
else:
    df = pd.read_csv('thermo_pyra_'+str(TC)+'_pseudo_eq_allData.csv')
        # List of thermo
for i in range(1,7,1):  
    cols = ['ΔE$_'+str(i)+'$', 'ΔE$_'+str(i)+'(ZPE)$', 'ΔG°$_'+str(i)+'$']
    v = df[cols].round(1)
    df[cols] = v
    df['combined_'+str(i)] = df[cols].apply(lambda row: '\n'.join(row.values.astype(str)), axis=1) 
    df.to_csv('df_with_values_to_copy.csv')   


END_SCRIPT

# Run script.py with python3
python3 gen_pseudo_rot_equil_csv_hist.py

cd ../
done

cd ../ # LEAVE IL FOLDER
done














