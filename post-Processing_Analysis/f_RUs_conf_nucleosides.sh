#!/bin/bash

#	f_RUs_conf_nucleosides.sh

#	09-jun-2023	Fixed titles of polar graphs in gen_summary_csv.py 	
# 	23-jan-2023	Created

# Description:
# This script obtains the conformation of the RU
# around the sugar at the glycosidic bond. It will classify
# the conformation as syn or anti.

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
## FUNCTION TO CALCULATE TORSION ANGLE GIVEN    ##
## 4 ATOMS (X, Y, Z) COORDINATES.               ##                                  
##                                              ##
##################################################

function calc_torsion_angle {

# Assign local variables to entries of the function
local a1_coord=$1
local a2_coord=$2
local a3_coord=$3
local a4_coord=$4

# Export the variables to be used in python
export a1_coord
export a2_coord
export a3_coord
export a4_coord

cat >calculate_torsion_value.py <<'END_SCRIPT'
# Import environmental variables from python libraries
import sys, os
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from numpy import loadtxt

# Import environmental variables from bash script
a1_coord = os.environ.get("a1_coord")
a2_coord = os.environ.get("a2_coord")
a3_coord = os.environ.get("a3_coord")
a4_coord = os.environ.get("a4_coord")

#print (a1_coord)
#print (a2_coord)

# Load each atom coordinates stored in files into numpy arrays variables
a1 = np.loadtxt(a1_coord, dtype=float)
a2 = np.loadtxt(a2_coord, dtype=float)
a3 = np.loadtxt(a3_coord, dtype=float)
a4 = np.loadtxt(a4_coord, dtype=float)

# Calculate the vector distances from a2 to a1, a3 to a2 and a4 to a3
v1 = a1 - a2
v2 = a2 - a3
v3 = a3 - a4

# Obtain the normalized vectors v1, v2 and v3
normalized_v1 = v1/np.linalg.norm(v1)
normalized_v2 = v2/np.linalg.norm(v2)
normalized_v3 = v3/np.linalg.norm(v3)

# Obtain the cross product v1 x v2 and v2 x v3 (normal vectors to the planes containing v1, v2 and v2, v3)
n1 = np.cross(normalized_v1, normalized_v2)
n2 = np.cross(normalized_v2, normalized_v3)

# Compute the sin of the torsion
m1 = np.cross(n1, normalized_v2)
sinTorsion = np.dot(m1, n2)
cosTorsion = np.dot(n1, n2)

#print ('a1', a1)
#print ('a2', a2)
#print ('v1', v1)
#print ('normalized v1', normalized_v1)

# Obtain the torsion angle
torsion_angle = np.degrees(np.arctan2(sinTorsion, cosTorsion))

# Convert angle to 0-360 scale if is negative
if (torsion_angle < 0):
    torsion_angle = 360 + torsion_angle

print (torsion_angle)

END_SCRIPT

# Run script.py with python3
torsion=$(python3 calculate_torsion_value.py $*)

}


########################################################
##                                                    ##
## FUNCTION TO OBTAINS THE BASE POSITION WITH RESPECT ##
## TO SUGAR AND CLASSIFY IN SYN OR ANTI.              ##
##                                                    ##
########################################################

function calc_RU_position {

# Load function input variables
local environment=$1
local TC=$2
local ring_conf=$3
local RU=$4

# Define input file for function that calculates the thermo equilibrium
log_file_prefix="$TC"_"$ring_conf"_"$RU"

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

# Cnvert log file to gjf to obtain XYZ coord
obabel $log_file -O cartesian_coordinates.gjf

# Delete initial lines, charge and multiplicity to leave only the cartesian coordinates and atoms
sed -i '1,6d' cartesian_coordinates.gjf

# Obtain only the x, y, z coordinates columns
awk '{print $2, $3, $4}' cartesian_coordinates.gjf | column -t > cartesian_coordinates.txt

################## Obtain the atoms that are involve in the X torsion angle #################
# Obtain the total number of atoms in the TC and the label for the hemiacetalic oxygen and the anomeric carbon
if [ "$TC" == "2deoxy_ribofuranose" ]; then
nAtoms=17
O_hemi=4
C_anomeric=8

elif [ "$TC" == "ribofuranose" ]; then 
nAtoms=18
O_hemi=4
C_anomeric=8

elif [ "$TC" == "2deoxy_ribopyranose" ]; then
nAtoms=17
O_hemi=6
C_anomeric=1

elif [ "$TC" == "ribopyranose" ]; then
nAtoms=18
O_hemi=6
C_anomeric=1

elif [ "$TC" == "threose" ]; then
nAtoms=14
O_hemi=2
C_anomeric=5
fi

# Obtain the label for the other two atoms of the RU
if [ "$RU" == "adenine" ]; then
a1_base=$((1+$nAtoms))
a2_base=$((9+$nAtoms))

elif [ "$RU" == "guanine" ]; then
a1_base=$((1+$nAtoms))
a2_base=$((11+$nAtoms))

elif [ "$RU" == "cytosine" ] || [ "$RU" == "thymine" ] || [ "$RU" == "uracil" ] || [ "$RU" == "CA" ] || [ "$RU" == "melamine" ]; then
a1_base=$((1+$nAtoms))
a2_base=$((2+$nAtoms))

elif [ "$RU" == "BA_cbond" ] || [ "$RU" == "BA_nbond" ] || [ "$RU" == "TARC_nbond" ]; then
a1_base=$((1+$nAtoms))
a2_base=$((2+$nAtoms))

elif [ "$RU" == "TARC_cbond" ]; then
a1_base=$((1+$nAtoms))
a2_base=$((5+$nAtoms))

fi

# Now copy lines with atoms in new file
sed -n ''$O_hemi'p' cartesian_coordinates.txt > a1_coord.txt
sed -n ''$C_anomeric'p' cartesian_coordinates.txt > a2_coord.txt
sed -n ''$a1_base'p' cartesian_coordinates.txt > a3_coord.txt
sed -n ''$a2_base'p' cartesian_coordinates.txt > a4_coord.txt

# Now pass this 4 files to function that calculates the torsion angle between the 4 atoms
calc_torsion_angle a1_coord.txt a2_coord.txt a3_coord.txt a4_coord.txt

# Print result in corresponding txt file
echo $torsion > torsionX_value.txt

export torsion

cat >obtain_base_position.py <<'END_SCRIPT'
# Import environmental variables from python libraries
import sys, os
import numpy as np 
import pandas as pd
from numpy import loadtxt

# Import environmental variables from bash script
torsion = os.environ.get("torsion")

# Ouput position of base with respect to sugar
summary = open("position_of_base.txt", "w")
if float(torsion) >= 90 and float(torsion) <= 270:    
    summary.write("Anti")
else:
    summary.write("Syn")
END_SCRIPT

# Run script.py with python3
python3 obtain_base_position.py

# Send results to general files outside in nucleosides folder
X_value=`cat torsionX_value.txt`
position=`cat position_of_base.txt`

echo $X_value >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/RU_conformation/$TC/chi_values.txt
echo $position >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/RU_conformation/$TC/TC_RU_relative_orientation.txt

echo $environment >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/RU_conformation/$TC/environment_list_names.txt
echo $TC >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/RU_conformation/$TC/TC_list_names.txt
echo $ring_conf >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/RU_conformation/$TC/ring_conf_list_names.txt
echo $RU >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/RU_conformation/$TC/RU_list_names.txt


}


##########################################################################################################################
############################################## WORK ON EACH FOLDER #######################################################
##########################################################################################################################

# Create folder for analysis and post-processing
mkdir $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/RU_conformation


# Open corresponding folder
cd nucleosides/final_gaussian_opt

for environment in "${environment[@]}"; do

cd $environment

for TC in "${TC_names[@]}"; do

cd $TC

if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

# Create corresponding folder for furanoses
mkdir $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/RU_conformation/$TC

for ring_conf_fura in "${Ring_Conf_furanose[@]}"; do

cd $ring_conf_fura

for RU in "${RU_names[@]}"; do

cd $RU

# Output message of status of script
echo "working on $environment/$TC/$ring_conf_fura/$RU"

# Run function to calculate the glycosidic torsion angle between RU and sugar 
calc_RU_position $environment $TC $ring_conf_fura $RU 




cd ../
done

cd ../
done

elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then

# Create corresponding folder for furanoses
mkdir $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/RU_conformation/$TC

for ring_conf_pyra in "${Ring_Conf_pyranose[@]}"; do

cd $ring_conf_pyra

for RU in "${RU_names[@]}"; do

cd $RU

# Output message of status of script
echo "working on $environment/$TC/$ring_conf_pyra/$RU"

# Run function to calculate the glycosidic torsion angle between RU and sugar 
calc_RU_position $environment $TC $ring_conf_pyra $RU 



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


####################################################################################################################################################
############################################## Generate final csv files with RUs conformations info ################################################
####################################################################################################################################################

cd post_proc_summary_results_csv_graphs/RU_conformation

cat >gen_summary_csv.py <<'END_SCRIPT'

# Import libraries
import sys, os
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from numpy import loadtxt
import matplotlib
import math
from scipy.stats import multivariate_normal
from scipy.stats import circvar
from scipy.stats import circmean
from scipy.stats import circstd
import scipy.stats as st
import matplotlib.patches as mpatches

environment = ["vacuum","water"]
TC_names = ["2deoxy_ribofuranose", "ribofuranose", "2deoxy_ribopyranose", "ribopyranose", "threose"]
alpha = ["2endo_alpha", "3endo_alpha", "D1C4_alpha", "D4C1_alpha"]
beta = ["2endo_beta", "3endo_beta", "D1C4_beta", "D4C1_beta"]

#####################################
#                                   #
# OBTAIN FREQUENCY TABLE FOR ANGLES #
#                                   #
#####################################

def frequencies_for_hist(angles_dataframe, lower_value, upper_value, step):
    # Create the bins or intervals
    bins = pd.interval_range(lower_value, upper_value, freq=step)    
    # Assign each value of angle to bin and count bin occurences, angles is the corresponding column of dataframe
    counts = pd.cut(angles_dataframe, bins).value_counts(sort=False)
    # Create frequency table
    freq_table = pd.Series(counts.values, index=bins)
    freq_table = pd.DataFrame({'Bins': freq_table.index, 'Count': freq_table.values})
    frequencies = np.array(freq_table['Count'])
    return (freq_table, frequencies)
    
    
    
# LOOP THROUGH THE DIFFERENT TCs FOLDERS

for TC in TC_names:

    os.chdir(TC) # open TC folders
    
    # Load all necesary txt files
    environment = np.loadtxt("environment_list_names.txt", dtype='str')
    TC = np.loadtxt("TC_list_names.txt", dtype='str')
    ring_conf = np.loadtxt("ring_conf_list_names.txt", dtype='str')
    RUs = np.loadtxt("RU_list_names.txt", dtype='str')
    position = np.loadtxt("TC_RU_relative_orientation.txt", dtype='str')
    chi_value = np.loadtxt("chi_values.txt", dtype='float')

    # Generate dataframe and csv file
    dataframe = pd.DataFrame({'Environment': environment,
         'TC': TC, 'Ring_conf': ring_conf, 'RU': RUs, 
         'Chi': chi_value, 'Position': position})
	 
    dataframe.to_csv("summary_db_RU_pos.csv", index=False)
    
    # Load the csv database
    db_all = pd.read_csv('summary_db_RU_pos.csv')
    
    # Obtain entries for vacuum and water environment from db
    db_vacuum = db_all[db_all['Environment'] == "vacuum"]
    db_water = db_all[db_all['Environment'] == "water"]
    
    # Obtain from db_env the entries for canonical and for not_canonical bases
    # List of RUs
    canonical_RU = ['adenine', 'guanine', 'cytosine', 'thymine', 'uracil']
    non_canonical_RU = ['TARC_cbond', 'TARC_nbond', 'BA_cbond', 'BA_nbond', 'melamine', 'CA']

    #vacuum
    db_vacuum_canonical = db_vacuum[db_vacuum['RU'].isin(canonical_RU)]
    chi_vacuum_canonical = db_vacuum_canonical['Chi']
    db_vacuum_non_canonical = db_vacuum[db_vacuum['RU'].isin(non_canonical_RU)]
    chi_vacuum_non_canonical = db_vacuum_non_canonical['Chi']

    # Water
    db_water_canonical = db_water[db_water['RU'].isin(canonical_RU)]
    chi_water_canonical = db_water_canonical['Chi']
    db_water_non_canonical = db_water[db_water['RU'].isin(non_canonical_RU)]
    chi_water_non_canonical = db_water_non_canonical['Chi']
    
    
    
    # Counts for alpha and beta
    # Vacuum
    all_vacuum_alpha = db_vacuum[db_vacuum['Ring_conf'].isin(alpha)]
    all_vacuum_beta = db_vacuum[db_vacuum['Ring_conf'].isin(beta)]   
    can_vacuum_alpha = db_vacuum_canonical[db_vacuum_canonical['Ring_conf'].isin(alpha)]
    can_vacuum_beta = db_vacuum_canonical[db_vacuum_canonical['Ring_conf'].isin(beta)]
    non_can_vacuum_alpha = db_vacuum_non_canonical[db_vacuum_non_canonical['Ring_conf'].isin(alpha)]
    non_can_vacuum_beta = db_vacuum_non_canonical[db_vacuum_non_canonical['Ring_conf'].isin(beta)]
    
    # Water
    all_water_alpha = db_water[db_water['Ring_conf'].isin(alpha)]
    all_water_beta = db_water[db_water['Ring_conf'].isin(beta)]   
    can_water_alpha = db_water_canonical[db_water_canonical['Ring_conf'].isin(alpha)]
    can_water_beta = db_water_canonical[db_water_canonical['Ring_conf'].isin(beta)]
    non_can_water_alpha = db_water_non_canonical[db_water_non_canonical['Ring_conf'].isin(alpha)]
    non_can_water_beta = db_water_non_canonical[db_water_non_canonical['Ring_conf'].isin(beta)]
    
    
    ################### Obtain values of phi2 angle #########################
    chi_all_vacuum_alpha = all_vacuum_alpha['Chi']
    chi_all_vacuum_beta = all_vacuum_beta['Chi']
    chi_can_vacuum_alpha = can_vacuum_alpha ['Chi']
    chi_can_vacuum_beta = can_vacuum_beta ['Chi']
    chi_non_can_vacuum_alpha = non_can_vacuum_alpha ['Chi']
    chi_non_can_vacuum_beta = non_can_vacuum_beta ['Chi']
    
    # Water
    chi_all_water_alpha = all_water_alpha['Chi']
    chi_all_water_beta = all_water_beta['Chi']
    chi_can_water_alpha = can_water_alpha ['Chi']
    chi_can_water_beta = can_water_beta ['Chi']
    chi_non_can_water_alpha = non_can_water_alpha ['Chi']
    chi_non_can_water_beta = non_can_water_beta ['Chi']
       
    ############### Obtain the frequency tables and counts ##################
    # Vacuum
    (freq_df_all_vacuum_alpha, count_chi_all_vacuum_alpha) = frequencies_for_hist(chi_all_vacuum_alpha, 
                                                                                  0, 360, 90)
    (freq_df_all_vacuum_beta, count_chi_all_vacuum_beta) = frequencies_for_hist(chi_all_vacuum_beta, 
                                                                                  0, 360, 90)  
                                                                                                                                                                                                                                                                                                                                           
    (freq_df_can_vacuum_alpha, count_chi_can_vacuum_alpha) = frequencies_for_hist(chi_can_vacuum_alpha, 
                                                                                  0, 360, 90)
    (freq_df_can_vacuum_beta, count_chi_can_vacuum_beta) = frequencies_for_hist(chi_can_vacuum_beta, 
                                                                                  0, 360, 90)                                                                             
    (freq_df_non_can_vacuum_alpha, count_chi_non_canonical_vacuum_alpha) = frequencies_for_hist(chi_non_can_vacuum_alpha,
                                                                                                0, 360, 90)
    (freq_df_non_can_vacuum_beta, count_chi_non_canonical_vacuum_beta) = frequencies_for_hist(chi_non_can_vacuum_beta,
                                                                                                0, 360, 90) 
                                                                                                                                                                                             
    # Water
    (freq_df_all_water_alpha, count_chi_all_water_alpha) = frequencies_for_hist(chi_all_water_alpha, 
                                                                                0, 360, 90)
    (freq_df_all_water_beta, count_chi_all_water_beta) = frequencies_for_hist(chi_all_water_beta, 
                                                                              0, 360, 90)  
                                                                                                                                                                                                                                                                                                                                           
    (freq_df_can_water_alpha, count_chi_can_water_alpha) = frequencies_for_hist(chi_can_water_alpha, 
                                                                                0, 360, 90)
    (freq_df_can_water_beta, count_chi_can_water_beta) = frequencies_for_hist(chi_can_water_beta, 
                                                                              0, 360, 90)                                                                             
    (freq_df_non_can_water_alpha, count_chi_non_canonical_water_alpha) = frequencies_for_hist(chi_non_can_water_alpha,
                                                                                              0, 360, 90)
    (freq_df_non_can_water_beta, count_chi_non_canonical_water_beta) = frequencies_for_hist(chi_non_can_water_beta,
                                                                                            0, 360, 90) 
                                                                                            
                                                                                            
                                                                                            
    # Obtain the circ mean, std and total datapoints 
    # Vacuum
    chi_can_vacuum_rad = np.radians(chi_vacuum_canonical)
    (circ_mean_vacuum_canonical, circ_std_vacuum_canonical) = (np.degrees(circmean(chi_can_vacuum_rad)),
                                                               np.degrees(circstd(chi_can_vacuum_rad)))
                                                                      
    chi_non_can_vacuum_rad = np.radians(chi_vacuum_non_canonical)    
    (circ_mean_vacuum_non_canonical, circ_std_vacuum_non_canonical) = (np.degrees(circmean(chi_non_can_vacuum_rad)),
                                                                       np.degrees(circstd(chi_non_can_vacuum_rad)))

    chi_all_vacuum_rad = np.radians(db_vacuum['Chi'])    
    (circ_mean_all_vacuum, circ_std_all_vacuum) = (np.degrees(circmean(chi_all_vacuum_rad)),
                                                                       np.degrees(circstd(chi_all_vacuum_rad)))

    # Water
    chi_can_water_rad = np.radians(chi_water_canonical)
    (circ_mean_water_canonical, circ_std_water_canonical) = (np.degrees(circmean(chi_can_water_rad)),
                                                             np.degrees(circstd(chi_can_water_rad)))
                                                                      
    chi_non_can_water_rad = np.radians(chi_water_non_canonical)    
    (circ_mean_water_non_canonical, circ_std_water_non_canonical) = (np.degrees(circmean(chi_non_can_water_rad)),
                                                                     np.degrees(circstd(chi_non_can_water_rad)))

    chi_all_water_rad = np.radians(db_water['Chi'])    
    (circ_mean_all_water, circ_std_all_water) = (np.degrees(circmean(chi_all_water_rad)),
                                                                     np.degrees(circstd(chi_all_water_rad)))
                                                                     
    ## Define vectors with all info for rose diagrams
    #Vacuum
    info_can_vacuum = [count_chi_can_vacuum_alpha, count_chi_can_vacuum_beta, circ_mean_vacuum_canonical, circ_std_vacuum_canonical]
    info_non_can_vacuum = [count_chi_non_canonical_vacuum_alpha, count_chi_non_canonical_vacuum_beta, circ_mean_vacuum_non_canonical, circ_std_vacuum_non_canonical]
    info_all_vacuum = [count_chi_all_vacuum_alpha, count_chi_all_vacuum_beta, circ_mean_all_vacuum, circ_std_all_vacuum]

    # Water
    info_can_water = [count_chi_can_water_alpha, count_chi_can_water_beta, circ_mean_water_canonical, circ_std_water_canonical]
    info_non_can_water = [count_chi_non_canonical_water_alpha, count_chi_non_canonical_water_beta, circ_mean_water_non_canonical, circ_std_water_non_canonical]
    info_all_water = [count_chi_all_water_alpha, count_chi_all_water_beta, circ_mean_all_water, circ_std_all_water]       
    
    
    
    ########## Create the rose diagrams #############
    fig, axn = plt.subplots(2, 3, sharey=False, sharex=False, subplot_kw=dict(projection="polar"), figsize=(20,14))
    for i, (ax, count_alpha, count_beta, mean, std) in enumerate([(axn.flat[0], info_can_vacuum[0], info_can_vacuum[1], 
                                                                   info_can_vacuum[2], info_can_vacuum[3]), 
                                                                  (axn.flat[1], info_non_can_vacuum[0], 
                                                                   info_non_can_vacuum[1], info_non_can_vacuum[2],
                                                                   info_non_can_vacuum[3]), 
                                                                  (axn.flat[2], info_all_vacuum[0], 
                                                                   info_all_vacuum[1], info_all_vacuum[2],
                                                                   info_all_vacuum[3]), 
                                                                  (axn.flat[3], info_can_water[0], info_can_water[1], 
                                                                   info_can_water[2], info_can_water[3]), 
                                                                  (axn.flat[4], info_non_can_water[0], 
                                                                   info_non_can_water[1], info_non_can_water[2],
                                                                   info_non_can_water[3]), 
                                                                  (axn.flat[5], info_all_water[0], info_all_water[1], 
                                                                   info_all_water[2], info_all_water[3])]):
        # Obtain the quadrants
        theta = np.radians(np.arange(0, 360, 90)) 
        mean = round(mean, 4)
        std = round(std, 4)
    
        # Obtain the range for freqs
        width = np.radians(90)  
    
        # Create bar plot
        bars = ax.bar(theta, count_alpha, width = width, edgecolor = 'black', align = 'edge', facecolor='gray', 
                      linewidth = 0.5, label = 'Alpha')
    
        bars2 = ax.bar(theta, count_beta, width = width, edgecolor = 'black', align = 'edge', facecolor='g', 
                       bottom = count_alpha, linewidth = 0.5, label = 'Beta')
    
        # Configure the grids
        ax.yaxis.grid(True,color='k',linestyle='--', linewidth=0.5)
        ax.xaxis.grid(True,color='black',linestyle='-', linewidth=0.8)    
        ax.xaxis.set_visible(True)
        ax.set_facecolor('xkcd:white')
        ax.grid(True)
    
        # Configure the axis ticks
        ax.tick_params(axis='x', colors='black', labelsize='22')
        ax.tick_params(axis='y', colors='blue', labelsize='22')
        plt.xticks(weight = 'bold')
        labelss = ax.get_yticklabels()
        for label in labelss:
            label.set_fontweight('bold')
    
        # Change position of y labels
        ax.set_rlabel_position(90)
    
        labels=['0°', "90°", '180°', "270°"]
        ax.set_xticks(theta)
        ax.set_xticklabels(labels, fontweight='bold') 
    
        x_tail = 0.0
        y_tail = 0.0
        x_head = float(mean)/180.*np.pi
        array = count_alpha + count_beta
        y_head = max(array)
        dx = x_head - x_tail
        dy = y_head - y_tail
        arrow = mpatches.FancyArrowPatch((x_tail, y_tail), (dx, dy),
                                          mutation_scale=30, facecolor = 'red', edgecolor = 'black')
        ax.add_patch(arrow)
    
        # Add the title to rose diagram
        if i == 0:
            ax.set_title('Canonical RUs\n\n(' + r'$\bar\chi$' + '=' + str(mean) + "°" +  "±" + str(std) + "°"
                         + ")" + '\n', fontweight="bold", fontsize=28)
        elif i == 1:
            ax.set_title('Non-Canonical RUs\n\n(' + r'$\bar\chi$' + '=' + str(mean) + "°" + "±" + str(std) + "°"
                         + ")" + '\n', fontweight="bold", fontsize=28)
        elif i == 2:
            ax.set_title('All RUs\n\n(' + r'$\bar\chi$' + '=' + str(mean) + "°" + "±" + str(std) + "°" + ")" + 
                         '\n', fontweight="bold", fontsize=28)
        elif i == 3:
            ax.set_title('\n(' + r'$\bar\chi$' + '=' + str(mean) + "°" + "±" + str(std) + "°" + ")" + 
                         '\n', fontweight="bold", fontsize=28)
        elif i == 4:
            ax.set_title('\n(' + r'$\bar\chi$' + '=' + str(mean) + "°" + "±" + str(std) + "°" + ")" 
                         + '\n', fontweight="bold", fontsize=28)
        elif i == 5:
            ax.set_title('\n(' + r'$\bar\chi$' + '=' + str(mean) + "°" + "±" + str(std) + "°" + ")" +
                         '\n', fontweight="bold", fontsize=28)
            ax.legend(['α', 'β'], fontsize=24, bbox_to_anchor=(1.3, 0.00), loc="lower right")
            leg = ax.get_legend()
            leg.legendHandles[0].set_color('grey')
            leg.legendHandles[1].set_color('green') 
     
        # little space between the bar and the label
        labelPadding = 4
        angles = np.radians(np.arange(45, 405, 90)) 
    
        # Set the bar numbers (First way)
        # Set the bar numbers (alpha anomers)
        for bar, angle, height, label in zip(bars, angles, count_alpha, count_alpha):
            # Labels are rotated. Rotation must be specified in degrees :(
            rotation = np.rad2deg(angle)        
            # Flip some labels upside down
            alignment = ""
            if angle >= np.pi/2 and angle < 3*np.pi/2:
                alignment = "right"
                rotation = rotation + 180
            else: 
                alignment = "left"
            if height != 0 and label != 0:
                ax.text(x=angle, 
                    y=bar.get_height(), 
                    s=label, 
                    ha=alignment, 
                    va='center', 
                    rotation=rotation, 
                    rotation_mode="anchor",
                    fontsize=18,
                    fontweight='bold')
                    
        # Set the bar numbers (First way)
        for bar, angle, height, label, ad_height in zip(bars2, angles, count_beta, count_beta, count_alpha):        
            # Labels are rotated. Rotation must be specified in degrees :(
            rotation = np.rad2deg(angle)        
            # Flip some labels upside down
            alignment = ""
            if angle >= np.pi/2 and angle < 3*np.pi/2:
                alignment = "right"
                rotation = rotation + 180
            else: 
                alignment = "left"
            if height != 0 and label != 0:
                ax.text(
                    x=angle, 
                    y=bar.get_height()+ad_height, 
                    s=label, 
                    ha=alignment, 
                    va='center', 
                    rotation=rotation, 
                    rotation_mode="anchor",
                    fontsize=18,
                    fontweight='bold')
       
    plt.tight_layout()       
    plt.savefig('Rose_diagram_'+str(TC[1])+'.pdf')                                                  
      
    os.chdir("..") # Leave TC folders

END_SCRIPT

# Run script.py with python3
python3 gen_summary_csv.py














