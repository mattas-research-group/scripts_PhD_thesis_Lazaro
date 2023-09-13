#!/bin/bash

# 12-Feb-2023  Created

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")

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
local log_file=$4

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
awk '{ printf "%.8e\n", $1*2625.5 }' zeropoint_energy_value.txt > zeropoint_energy_kj.txt

# Make copy of this file with appropiate name for component
cp zeropoint_energy_kj.txt "$log_file"_zeropoint_energy_kj.txt

# Assign variable to zero point energy in kcal/mol txt file
zero_point_energy=`cat "$log_file"_zeropoint_energy_kj.txt`

############################################# OBTAIN CORRECTED ENERGY ################################
line_corrected_energy=$(($line_n_zeropoint_energy+1))
sed -n ''$line_corrected_energy'p' "$log_file" > corrected_energy_react.txt

# Print 2nd column with actual energy value
awk '{ print $7 }' corrected_energy_react.txt > corrected_energy_value_react.txt

# Convert energy to kcal/mol
awk '{ printf "%.8e\n", $1*2625.5 }' corrected_energy_value_react.txt > corrected_energy_kj.txt

# Make copy of this file with appropiate name for component
cp corrected_energy_kj.txt "$log_file"_corrected_energy_kj.txt

# Assign variable to zero point energy in kcal/mol txt file
corrected_energy=`cat "$log_file"_corrected_energy_kj.txt`

############################################# OBTAIN ENTHALPY ################################
line_enthalpy=$(($line_n_zeropoint_energy+2))
sed -n ''$line_enthalpy'p' "$log_file" > entalphy_react.txt

# Print 2nd column with actual energy value
awk '{ print $7 }' entalphy_react.txt > enthalpy_value_react.txt

# Convert energy to kcal/mol
awk '{ printf "%.8e\n", $1*2625.5 }' enthalpy_value_react.txt > entalphy_kj.txt

# Make copy of this file with appropiate name for component
cp entalphy_kj.txt "$log_file"_entalphy_kj.txt

# Assign variable to zero point energy in kcal/mol txt file
enthalpy=`cat "$log_file"_entalphy_kj.txt`

############################################ OBTAIN THE FREE ENERGY ##################################
line_dG=$(($line_n_zeropoint_energy+3))
sed -n ''$line_dG'p' "$log_file" > dG_react.txt

# Print 2nd column with actual energy value
awk '{ print $8 }' dG_react.txt > dG_value_react.txt

# Convert energy to kcal/mol
awk '{ printf "%.8e\n", $1*2625.5 }' dG_value_react.txt > G_kj.txt

# Make copy of this file with appropiate name for component
cp G_kj.txt "$log_file"_G_kj.txt

# Assign variable to zero point energy in kcal/mol txt file
free_energy=`cat "$log_file"_G_kj.txt`

# OUTPUT ENERGIES AND OTHER INFO IN CORRESPONDING SUMMARY FILES IN ANALYSIS FOLDER
if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

echo $environment >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/furanoses/"$ring_conf"_environment_list.txt
echo $TC >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/furanoses/"$ring_conf"_TCs_names_list.txt
echo $ring_conf >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/furanoses/"$ring_conf"_ring_conf_names_list.txt

echo $zero_point_energy >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/furanoses/"$ring_conf"_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/furanoses/"$ring_conf"_corrected_energies.txt
echo $enthalpy >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/furanoses/"$ring_conf"_enthalpies.txt
echo $free_energy >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/furanoses/"$ring_conf"_free_energies.txt

else

echo $environment >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/pyranoses/"$ring_conf"_environment_list.txt
echo $TC >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/pyranoses/"$ring_conf"_TCs_names_list.txt
echo $ring_conf >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/pyranoses/"$ring_conf"_ring_conf_names_list.txt

echo $zero_point_energy >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/pyranoses/"$ring_conf"_zero_point_energies.txt
echo $corrected_energy >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/pyranoses/"$ring_conf"_corrected_energies.txt
echo $enthalpy >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/pyranoses/"$ring_conf"_enthalpies.txt
echo $free_energy >> $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/pyranoses/"$ring_conf"_free_energies.txt

fi

}

##########################################################################################################################
############################################## WORK ON EACH FOLDER #######################################################
##########################################################################################################################

# Create folder for analysis and post-processing
mkdir $MD_RUN/TC/post_proc_summary_results_csv_graphs
mkdir $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium

# Open corresponding folder
cd TC
cd final_gaussian_opt

for environment in "${environment[@]}"; do

cd $environment

for TC in "${TC_names[@]}"; do

cd $TC

if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

mkdir $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/furanoses

for ring_conf_fura in "${Ring_Conf_furanose[@]}"; do

cd $ring_conf_fura

echo "working on "$environment"/"$TC"/"$ring_conf_fura""

# Run function to extract thermo parameters from log files for nucleoside, RU and TC 
extract_thermo_params_pseu_rot $environment $TC $ring_conf_fura "$TC"_"$ring_conf_fura".log



cd ../ # LEAVE CURRENT RING CONF

done

else

mkdir $MD_RUN/TC/post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium/pyranoses

for ring_conf_pyra in "${Ring_Conf_pyranose[@]}"; do

cd $ring_conf_pyra

echo "working on "$environment"/"$TC"/"$ring_conf_pyra""

# Run function to extract thermo parameters from log files for nucleoside, RU and TC 
extract_thermo_params_pseu_rot $environment $TC $ring_conf_pyra "$TC"_"$ring_conf_pyra".log




cd ../ # LEAVE CURRENT RING CONF

done

fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT

cd ../ # LEAVE THE FINAL_GAUSSIAN_OPT FOLDER

###########################################################################################################################################
############################################## WORK ON CONDENSATION REACTION FOLDER #######################################################
###########################################################################################################################################

cd post_proc_summary_results_csv_graphs/sugar_ring_pseudorot_equilibrium

# Create and run python script to gen summary csv files and freq hist 
cat >gen_pseudo_rot_equil_csv_hist.py <<'END_SCRIPT'
# Import environmental variables from python libraries
import sys, os
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from numpy import loadtxt
import seaborn as sb
import scipy as scipy
from scipy import stats

# Create necesary stack of options
environment = ['vacuum','water']
TC_names = ["2deoxy_ribofuranose", "ribofuranose", "2deoxy_ribopyranose", "ribopyranose", "threose"]
Ring_Conf_furanose = ["2endo_alpha", "2endo_beta", "3endo_alpha", "3endo_beta"]
Ring_Conf_pyranose = ["D1C4_alpha", "D1C4_beta", "D4C1_alpha", "D4C1_beta"]
thermo_list = ["zero_point_energies", "corrected_energies", "enthalpies", "free_energies"]
thermo_params = ['ΔzpE', 'ΔcorrE', 'ΔH', 'ΔG']
sugars = ["furanoses", "pyranoses"]


def obtain_final_db (db, thermo): 
    dG_list = []
    for i in range(1, 7, 1):
        G = str(thermo)+"#"+str(i)
        dG_list.append(G)
    db = db.filter(dG_list, axis=1)
    return db    






for ring_size in sugars:

    os.chdir(ring_size)
    
    if ring_size == "furanoses":
    
        # Create initial dataframe
        env = np.loadtxt("2endo_alpha_environment_list.txt", dtype='str')
        TC = np.loadtxt("2endo_alpha_TCs_names_list.txt", dtype='str')
        dataframe = pd.DataFrame({'environment': env, 'TC': TC})   
        
        # Calculate the dS parameters
        for ring_conf in Ring_Conf_furanose:
    
            for i, params in enumerate(thermo_list):
        
                # Get label for param in the final dataframe
                thermo = thermo_params[i]         
                # Output given param in the final dataframe
                param = np.loadtxt(str(ring_conf)+'_'+str(params)+'.txt', dtype='float')
                dataframe[str(ring_conf)+'_'+str(thermo)] = param 
                dataframe[str(ring_conf)+'_TΔS'] = dataframe[str(ring_conf) + '_ΔH'] - dataframe[str(ring_conf) + '_ΔG']
            
        # Obtain the ddparams
        therm_params = ['ΔzpE', 'ΔcorrE', 'ΔH', 'ΔG', 'TΔS']
        for i, params in enumerate(therm_params):
    
            dataframe[str(params)+'#1'] = dataframe['2endo_alpha_'+str(params)] - dataframe['3endo_alpha_'+str(params)]
            dataframe[str(params)+'#2'] = dataframe['2endo_beta_'+str(params)] - dataframe['2endo_alpha_'+str(params)]
            dataframe[str(params)+'#3'] = dataframe['2endo_beta_'+str(params)] - dataframe['3endo_beta_'+str(params)]
            dataframe[str(params)+'#4'] = dataframe['3endo_beta_'+str(params)] - dataframe['3endo_alpha_'+str(params)]        
            dataframe[str(params)+'#5'] = dataframe['2endo_beta_'+str(params)] - dataframe['3endo_alpha_'+str(params)]        
            dataframe[str(params)+'#6'] = dataframe['3endo_beta_'+str(params)] - dataframe['2endo_alpha_'+str(params)]   
    
            # Select the vacuum and water rows
            db_vacuum = dataframe[dataframe['environment'] == "vacuum"]
            db_water = dataframe[dataframe['environment'] == "water"]
            
            # Change names
	    TC_old_names = ["2deoxy_ribofuranose", "ribofuranose", "threose"]
            TC_new_names = ["2dRibf", "Ribf", "Tho"]
            
            for i, TC_n in enumerate(TC_old_names):
                db_water['TC'] = db_water['TC'].str.replace(TC_n, TC_new_names[i])
            
            # Generate from this df the columns for params values + RU + environment
            vacuum_2deoxy_ribofuranose = obtain_final_db (db_vacuum[db_vacuum['TC'] == "2deoxy_ribofuranose"], params)
            vacuum_ribofuranose = obtain_final_db (db_vacuum[db_vacuum['TC'] == "ribofuranose"], params)
            vacuum_threose = obtain_final_db (db_vacuum[db_vacuum['TC'] == "threose"], params)

            water_2deoxy_ribofuranose = obtain_final_db (db_water[db_water['TC'] == "2dRibf"], params)
            water_ribofuranose = obtain_final_db (db_water[db_water['TC'] == "Ribf"], params)
            water_threose = obtain_final_db (db_water[db_water['TC'] == "Tho"], params)
            
            # Generate the bar plot
            fig, axn = plt.subplots(2, 3, sharex=False, sharey=False, figsize=(16,12)) 
            for i, (ax, dataframe, env, param) in enumerate([(axn.flat[0], dG_vacuum_2deoxy_ribofuranose, "vacuum", params), 
                                                           (axn.flat[1], dG_vacuum_ribofuranose, "vacuum", params), 
                                                           (axn.flat[2], dG_vacuum_threose, "vacuum", params),
                                                           (axn.flat[3], dG_water_2deoxy_ribofuranose, "water", params),
                                                           (axn.flat[4], dG_water_ribofuranose, "water", params),
                                                           (axn.flat[5], dG_water_threose, "water", params)
                                                           ]):
                # Create plot
                bar = dataframe.plot(kind='bar', ax=ax, colormap='Greys', width=0.8, edgecolor = "black", xlabel='', 
                                     legend=False, fontsize=20)
                # Delete x ticks
                bar.set_xticks([]) 
    
                # Create horizontal ax line
                bar.axhline(y = 0.0, color = 'black', linestyle = '--')
    
                # Add x axis to bar plots
                if ax == axn.flat[3]:
                    bar.set_xlabel("2deoxy_ribofuranose", fontsize=16, fontdict=dict(weight='bold'))
                elif ax == axn.flat[4]:
                    bar.set_xlabel("ribofuranose", fontsize=16, fontdict=dict(weight='bold'))
                elif ax == axn.flat[5]:
                    bar.set_xlabel("threose", fontsize=16, fontdict=dict(weight='bold'))
                    bar.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=16)
        
                # Add y axis to bar plots
                if ax == axn.flat[0] or ax == axn.flat[3]:
                    bar.set_ylabel(str(param)+"$^\circ$(kJ/mol)", fontsize=16, fontdict=dict(weight='bold'))
        
                # Add titles
                if ax == axn.flat[0] or ax == axn.flat[1] or ax == axn.flat[2]:
                    bar.set_title('Vacuum', weight='bold', fontsize=18)
                elif ax == axn.flat[3] or ax == axn.flat[4] or ax == axn.flat[5]:
                    bar.set_title('Water', weight='bold', fontsize=18)
        
                # Add containers for the bars
                for container in bar.containers:
                    bar.bar_label(container, fmt='%.1f', fontsize=16)
               
            plt.tight_layout()    
            plt.savefig('pseudorot_equil_furanoses_' + str(params) + '_graph1.pdf')
            
            #Obtain the dataframes for each thermo param in vacuum and water
            dG_vacuum = obtain_final_db (db_vacuum)
            dG_water = obtain_final_db (db_water)
            
            # Generate graph
            fig, axn = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(8,10)) 
            for i, (ax, dataframe, env) in enumerate([(axn.flat[0], dG_vacuum, "vacuum"), 
                                                      (axn.flat[1], dG_water, "water"), 
                                                    ]):
    
                bar = dataframe.plot(kind='bar', ax=ax, colormap='Greys', rot = 360, width=0.8, edgecolor = "black", xlabel='', 
                                     legend=False, fontsize=20)
    
                bar.axhline(y = 0.0, color = 'black', linestyle = '--')
        
                # Add y axis to bar plot
                bar.set_ylabel("ΔG$^\circ$(kJ/mol)", fontsize=16, fontdict=dict(weight='bold'))
    
                # Set x ticks for bars and legend
                if ax == axn.flat[1]:
                    bar.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=16)
                    labels = bar.get_xticklabels() 
                    [label.set_fontweight('bold') for label in labels]
                else:
                    bar.set_xticks([]) 
    
                # Add titles
                if ax == axn.flat[0]:
                    bar.set_title('Vacuum', weight='bold', fontsize=18)
                elif ax == axn.flat[1]:
                    bar.set_title('Water', weight='bold', fontsize=18)
        
                # Add containers for the bars
                for container in bar.containers:
                bar.bar_label(container, fmt='%.1f', fontsize=16)
                   
            plt.tight_layout()    
            plt.savefig('pseudorot_equil_furanoses_' + str(params) + '_graph2.pdf')
            
            
    else:    # Work on pyranoses
    
        # Create initial dataframe
        env = np.loadtxt("D1C4_alpha_environment_list.txt", dtype='str')
        TC = np.loadtxt("D1C4_alpha_TCs_names_list.txt", dtype='str')
        dataframe = pd.DataFrame({'environment': env, 'TC': TC})   
        
        # Calculate the dS parameters
        for ring_conf in Ring_Conf_pyranose:
    
            for i, params in enumerate(thermo_list):
        
                # Get label for param in the final dataframe
                thermo = thermo_params[i]         
                # Output given param in the final dataframe
                param = np.loadtxt(str(ring_conf)+'_'+str(params)+'.txt', dtype='float')
                dataframe[str(ring_conf)+'_'+str(thermo)] = param 
                dataframe[str(ring_conf)+'_TΔS'] = dataframe[str(ring_conf) + '_ΔH'] - dataframe[str(ring_conf) + '_ΔG']
            
        # Obtain the ddparams
        therm_params = ['ΔzpE', 'ΔcorrE', 'ΔH', 'ΔG', 'TΔS']
        for i, params in enumerate(therm_params):
    
            dataframe[str(params)+'#1'] = dataframe['D4C1_alpha_'+str(params)] - dataframe['D1C4_alpha_'+str(params)] ##### CONTINUE HEREEEEEEEEEEEEEEEEEE !!!!!!!!!!!!!!!!
            dataframe[str(params)+'#2'] = dataframe['D4C1_beta_'+str(params)] - dataframe['D4C1_alpha_'+str(params)]
            dataframe[str(params)+'#3'] = dataframe['D4C1_beta_'+str(params)] - dataframe['D1C4_beta_'+str(params)]
            dataframe[str(params)+'#4'] = dataframe['D1C4_beta_'+str(params)] - dataframe['D1C4_alpha_'+str(params)]        
            dataframe[str(params)+'#5'] = dataframe['D4C1_beta_'+str(params)] - dataframe['D1C4_alpha_'+str(params)]        
            dataframe[str(params)+'#6'] = dataframe['D4C1_alpha_'+str(params)] - dataframe['D1C4_beta_'+str(params)]   
    
            # Select the vacuum and water rows
            db_vacuum = dataframe[dataframe['environment'] == "vacuum"]
            db_water = dataframe[dataframe['environment'] == "water"]
            
            # Change names
	    TC_old_names = ["2deoxy_ribopyranose", "ribopyranose"]
            TC_new_names = ["2dRib", "Rib"]
            
            for i, TC_n in enumerate(TC_old_names):
                db_water['TC'] = db_water['TC'].str.replace(TC_n, TC_new_names[i])
            
            # Generate from this df the columns for params values + RU + environment
            vacuum_2deoxy_ribopyranose = obtain_final_db (db_vacuum[db_vacuum['TC'] == "2dRib"], params)
            vacuum_ribopyranose = obtain_final_db (db_vacuum[db_vacuum['TC'] == "Rib"], params)

            water_2deoxy_ribopyranose = obtain_final_db (db_water[db_water['TC'] == "2dRib"], params)
            water_ribopyranose = obtain_final_db (db_water[db_water['TC'] == "Rib"], params)
            
            # Generate the bar plot
            fig, axn = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(16,12)) 
            for i, (ax, dataframe, env, param) in enumerate([(axn.flat[0], dG_vacuum_2deoxy_ribopyranose, "vacuum", params), 
                                                           (axn.flat[1], dG_vacuum_ribopyranose, "vacuum", params), 
                                                           (axn.flat[2], dG_water_2deoxy_ribopyranose, "water", params),
                                                           (axn.flat[3], dG_water_ribopyranose, "water", params)
                                                           ]):
                # Create plot
                bar = dataframe.plot(kind='bar', ax=ax, colormap='Greys', width=0.8, edgecolor = "black", xlabel='', 
                                     legend=False, fontsize=20)
                # Delete x ticks
                bar.set_xticks([]) 
    
                # Create horizontal ax line
                bar.axhline(y = 0.0, color = 'black', linestyle = '--')
    
                # Add x axis to bar plots
                if ax == axn.flat[2]:
                    bar.set_xlabel("2deoxy_ribopyranose", fontsize=16, fontdict=dict(weight='bold'))
                elif ax == axn.flat[3]:
                    bar.set_xlabel("ribopyranose", fontsize=16, fontdict=dict(weight='bold'))
                    bar.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=16)
        
                # Add y axis to bar plots
                if ax == axn.flat[0] or ax == axn.flat[2]:
                    bar.set_ylabel(str(param)+"$^\circ$(kJ/mol)", fontsize=16, fontdict=dict(weight='bold'))
        
                # Add titles
                if ax == axn.flat[0] or ax == axn.flat[1]:
                    bar.set_title('Vacuum', weight='bold', fontsize=18)
                elif ax == axn.flat[2] or ax == axn.flat[3]:
                    bar.set_title('Water', weight='bold', fontsize=18)
        
                # Add containers for the bars
                for container in bar.containers:
                    bar.bar_label(container, fmt='%.1f', fontsize=16)
               
            plt.tight_layout()    
            plt.savefig('pseudorot_equil_pyranoses_' + str(params) + '_graph1.pdf') 
            
            #Obtain the dataframes for each thermo param in vacuum and water
            dG_vacuum = obtain_final_db (db_vacuum)
            dG_water = obtain_final_db (db_water)
            
            # Generate graph
            fig, axn = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(8,10)) 
            for i, (ax, dataframe, env) in enumerate([(axn.flat[0], dG_vacuum, "vacuum"), 
                                                      (axn.flat[1], dG_water, "water"), 
                                                    ]):
    
                bar = dataframe.plot(kind='bar', ax=ax, colormap='Greys', rot = 360, width=0.8, edgecolor = "black", xlabel='', 
                                     legend=False, fontsize=20)
    
                bar.axhline(y = 0.0, color = 'black', linestyle = '--')
        
                # Add y axis to bar plot
                bar.set_ylabel("ΔG$^\circ$(kJ/mol)", fontsize=16, fontdict=dict(weight='bold'))
    
                # Set x ticks for bars and legend
                if ax == axn.flat[1]:
                    bar.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=16)
                    labels = bar.get_xticklabels() 
                    [label.set_fontweight('bold') for label in labels]
                else:
                    bar.set_xticks([]) 
    
                # Add titles
                if ax == axn.flat[0]:
                    bar.set_title('Vacuum', weight='bold', fontsize=18)
                elif ax == axn.flat[1]:
                    bar.set_title('Water', weight='bold', fontsize=18)
        
                # Add containers for the bars
                for container in bar.containers:
                bar.bar_label(container, fmt='%.1f', fontsize=16)
                   
            plt.tight_layout()    
            plt.savefig('pseudorot_equil_pyranoses_' + str(params) + '_graph2.pdf')
    
     
    os.chdir("..") # Leave ring_type folders



END_SCRIPT

# Run script.py with python3
python3 gen_pseudo_rot_equil_csv_hist.py













