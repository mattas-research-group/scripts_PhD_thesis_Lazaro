#!/bin/bash

#	f_CPpuckering_nucleosides.sh

#	09-jun-2023	Fixed titles of polar graphs in gen_summary_csv.py 	
# 	23-jan-2023	Created

# Description:
# This script calculates the Cremer-Pople parameters for 
# nucleosides.

# Note on modifying CP kernel surfaces for pyranoses and circular plots for fura:
# Lines to change to edit polar plots:
#line #957
#
# Lines to change to edit kernel surfaces:
# line 195 and 196: to convert scale of phi from -180--180 to 0--360
# 
# line #1211 to assign new scale to phi in surf 
# line #1356-1361 to invert y axis and adjust labels for 1C4 and 4C1 

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a RU_names=("adenine" "guanine" "cytosine" "thymine" "uracil" "TARC_cbond" "TARC_nbond" "BA_cbond" "BA_nbond" "CA" "melamine")



###################################################################
##                                                               ##
## Function to obtain the generalized puckering parameters from  ##
## Cremer & Pople  (Cremer&Pople (1975)) for 5- and 6-member     ##
## rings.                                                        ##                                                 
##                                                               ##
###################################################################

function calc_puckering_param {

# Load function input variables
local environment=$1
local TC=$2
local ring_conf=$3
local RU=$4

export TC

cat >calc_ring_puckering.py <<'END_SCRIPT'

# Import environmental variables from python libraries
import sys, os
import numpy as np 
import pandas as pd
from numpy import loadtxt

# Import environmental variables from bash script for angles in degrees
TC = os.environ.get("TC")

###### Obtain the Cremer & Pople parameters #######

# Necesary functions
################################
########## Fix Zero ############
################################

def fixzero(x):
    x_ = np.array([0.0]) if np.allclose(0,x, rtol=1e-06, atol=1e-08) else x 
    return x_ 

def GetMeanPlane(coordinates):
    """
    Compute the mean plane
    
    Input:
    
    coordinates: array
    Return:
    
    R1: array
    R2: array
    """
    N = coordinates.shape[0] # ring size
    R1 = np.dot(np.sin(2*np.pi*np.arange(0,N)/N),coordinates)
    R2 = np.dot(np.cos(2*np.pi*np.arange(0,N)/N),coordinates)
    return R1, R2
    
def GetNormal(coordinates):
    """
    Compute normal to mean plan
    Input:
    coordinates: array
    Output:
    unit_normal: array
    """
    R1,R2 = GetMeanPlane(coordinates)
    cross_product = np.cross(R1,R2)
    unit_normal = cross_product/np.linalg.norm(cross_product) 
    return unit_normal
    
def Displacement(coordinates):
    """
    Compute the displacement (z) 
    
    Input:
    
    coordinates: array
    Output:
    Z: array
    """
    n = GetNormal(coordinates)
    z = np.dot(coordinates, n)
    return z
    
def GetRingPuckerCoords(coordinates):

    """
    Compute Ring Pucker Parameters
    
    Input:
    coordinates: array
    
    Return:
    qs: puckering amplitude (q_i>=0 for all i)
    
    angle: angle defined in 0<= phi_i <= 2pi 
    """

    N = coordinates.shape[0]  # number of atoms in the ring
    z = Displacement(coordinates)
    if N>4 and N<=20: # In our analysis, we fit it to be smaller than 16.
        if (N%2 == 0): # N even
            print('this has even atoms ring, e.g. 6')
            m = range(2,int((N/2)))
            cos_component = [np.dot(z,np.cos(2*np.pi*k*np.arange(0,N)/N)) for k in m]
            sin_component = [np.dot(z,np.sin(2*np.pi*k*np.arange(0,N)/N)) for k in m]
            qcos = fixzero(np.sqrt(2/N)*np.array(cos_component))
            qsin = fixzero(-np.sqrt(2/N)*np.array(sin_component))
            q = np.sqrt(qsin**2 + qcos**2)
            amplitude = np.append(q, (1/np.sqrt(N))*np.dot(z,np.cos(np.arange(0,N)*np.pi)).sum()).tolist()
            phase_angle = np.arctan2(qsin,qcos).tolist()
            Q = np.sqrt(np.square(amplitude).sum())
        else: # N odd
            print('this has odd atoms ring, e.g. 5')
            m = range(2,int((N-1)/2)+1)
            cos_component = [np.dot(z,np.cos(2*np.pi*k*np.arange(0,N)/N)) for k in m]
            sin_component = [np.dot(z,np.sin(2*np.pi*k*np.arange(0,N)/N)) for k in m]
            qcos = fixzero(np.sqrt(2/N)*np.array(cos_component))
            qsin = fixzero(-np.sqrt(2/N)*np.array(sin_component))
            amplitude = np.sqrt(qsin**2 + qcos**2).tolist()
            phase_angle = np.arctan2(qsin,qcos).tolist()            
            Q = np.sqrt(np.square(amplitude).sum())
    else:
        print("Ring Size is too big or too small")
    return amplitude, Q, phase_angle, N

# Loop through the different TCs to obtain necesary info
if (TC == '2deoxy_ribofuranose') or (TC == 'ribofuranose') or (TC == 'threose'):

    # Load each atom coordinates into numpy array
    R1 = np.loadtxt("O4_fura_ring.txt", dtype=float)
    R2 = np.loadtxt("C1_fura_ring.txt", dtype=float)
    R3 = np.loadtxt("C2_fura_ring.txt", dtype=float)
    R4 = np.loadtxt("C3_fura_ring.txt", dtype=float)
    R5 = np.loadtxt("C4_fura_ring.txt", dtype=float)
    
    # Obtain general matrix with coordinates
    coordinates_matrix = np.array([R1, R2, R3, R4, R5])

elif (TC == '2deoxy_ribopyranose') or (TC == 'ribopyranose'):

    # Load each atom coordinates into numpy array
    R1 = np.loadtxt("O4_pyra_ring.txt", dtype=float)
    R2 = np.loadtxt("C1_pyra_ring.txt", dtype=float)
    R3 = np.loadtxt("C2_pyra_ring.txt", dtype=float)
    R4 = np.loadtxt("C3_pyra_ring.txt", dtype=float)
    R5 = np.loadtxt("C4_pyra_ring.txt", dtype=float)
    R6 = np.loadtxt("C5_pyra_ring.txt", dtype=float)

    # Obtain general matrix with coordinates
    coordinates_matrix = np.array([R1, R2, R3, R4, R5, R6])      
  
# Run functions to generate ring puckering
[amplitude, Q, phase_angle, N] = GetRingPuckerCoords(coordinates_matrix)

if (N == 6): ########## CASE OF PYRANOSE RINGS

    theta = np.degrees(np.arccos(amplitude[1]/Q))
    #if theta < 0:
        #theta = 360 + theta
        
    phase_angle = np.degrees(phase_angle)
    #if phase_angle < 0:
        #phase_angle = 360 + phase_angle
        
    ################################################# 
    ## Ouput all parameters in corresponding files ##
    #################################################    

    # Ouput tm and P in different files
    summary_amplitude_q2 = open("amplitude_q2_from_C_P.txt", "w")
    summary_amplitude_q3 = open("amplitude_q3_from_C_P.txt", "w")
    summary_total_amplitude = open("total_amplitude_Q_from_C_P.txt", "w")
    summary_phase_angle = open("phase_angle_Phi2.txt", "w")
    summary_theta_angle = open("phase_angle_theta.txt", "w")
    summary_amplitude_q2.write(str(amplitude[0]))
    summary_amplitude_q3.write(str(amplitude[1]))
    summary_total_amplitude.write(str(Q))
    summary_phase_angle.write(str(phase_angle))
    summary_theta_angle.write(str(theta))

    ###################
    ## Classify ring ##
    ###################

    # Create conditions to output 6-member ring conformation base on phase angle value (Theta)
    summary_pyra_ring_conf_from_C_P = open("6-member_ring_conformation_C_P_from_theta.txt", "w")
    
    # Classify 6-member ring

    if theta == 0:
        summary_pyra_ring_conf_from_C_P.write("4C1_chair")
    elif theta > 0 and theta < 50.8:
        summary_pyra_ring_conf_from_C_P.write("4C1_chair---half-chair")
    elif theta == 50.8:
        summary_pyra_ring_conf_from_C_P.write("Half-chair")
    elif theta > 50.8 and theta < 54.7:
        summary_pyra_ring_conf_from_C_P.write("half-chair---envelope")
    elif theta == 54.7:
        summary_pyra_ring_conf_from_C_P.write("Envelope")
    elif theta > 54.7 and theta < 90:
        summary_pyra_ring_conf_from_C_P.write("envelope---boat")
    elif theta == 90:
        summary_pyra_ring_conf_from_C_P.write("Boat")
    elif theta > 90 and theta < 125.3:
        summary_pyra_ring_conf_from_C_P.write("boat---envelope")
    elif theta == 125.3:
        summary_pyra_ring_conf_from_C_P.write("Envelope")
    elif theta > 125.4 and theta < 129.2:
        summary_pyra_ring_conf_from_C_P.write("envelope---half-chair")
    elif theta == 129.2:
        summary_pyra_ring_conf_from_C_P.write("Half-chair")
    elif theta > 129.2 and theta < 180:
        summary_pyra_ring_conf_from_C_P.write("half-chair---1C4_chair")
    elif theta == 180:
        summary_pyra_ring_conf_from_C_P.write("1C4_chair")
	
elif (N == 5): ########## CASE OF FURANOSE RINGS

    phase_angle = np.degrees(phase_angle)
    if phase_angle < 0:
        phase_angle = 360 + phase_angle

    ################################################# 
    ## Ouput all parameters in corresponding files ##
    #################################################

    # Ouput tm and P in different files
    summary_amplitude = open("amplitude_q2_from_C_P.txt", "w")
    summary_total_amplitude = open("total_amplitude_Q_from_C_P.txt", "w")
    summary_phase_angle = open("phase_angle_Phi2.txt", "w")
    summary_amplitude.write(str(amplitude))
    summary_total_amplitude.write(str(Q))
    summary_phase_angle.write(str(phase_angle))

    ###################
    ## Classify ring ##
    ###################

    # Create conditions to output 5-member ring conformation base on phase angle (Phi2) value
    summary_fura_ring_conf_from_C_P = open("5-member_ring_conformation_C_P.txt", "w")

    # Classify 5-member ring accordingly to Altona & Sundaralingam phase angle (P)

    if phase_angle == 0:
        summary_fura_ring_conf_from_C_P.write("0E(O4-endo)")
    elif phase_angle > 0 and phase_angle < 18:
        summary_fura_ring_conf_from_C_P.write("0E---0T1")
    elif phase_angle == 18:
        summary_fura_ring_conf_from_C_P.swrite("0T1")
    elif phase_angle > 18 and phase_angle < 36:
        summary_fura_ring_conf_from_C_P.write("0T1---E1")
    elif phase_angle == 36:
        summary_fura_ring_conf_from_C_P.swrite("E1(C1'-exo)")
    elif phase_angle > 36 and phase_angle < 54:
        summary_fura_ring_conf_from_C_P.write("E1---2T1")
    elif phase_angle == 54:
        summary_fura_ring_conf_from_C_P.swrite("2T1")    
    elif phase_angle > 54 and phase_angle < 72:
        summary_fura_ring_conf_from_C_P.write("2T1--2E")
    elif phase_angle == 72:
        summary_fura_ring_conf_from_C_P.swrite("2E(C2'-endo)(South)")        
    elif phase_angle > 72 and phase_angle < 90:
        summary_fura_ring_conf_from_C_P.write("2E---2T3(South)")
    elif phase_angle == 90:
        summary_fura_ring_conf_from_C_P.swrite("2T3(C2'-endo_C3'-exo)(South)")    
    elif phase_angle > 90 and phase_angle < 108:
        summary_fura_ring_conf_from_C_P.write("2T3---E3(South)")
    elif phase_angle == 108:
        summary_fura_ring_conf_from_C_P.write("E3(C3'-exo)(South)")    
    elif phase_angle > 108 and phase_angle < 126:
        summary_fura_ring_conf_from_C_P.write("E3---4T3(C1'-exo)")
    elif phase_angle == 126:
        summary_fura_ring_conf_from_C_P.write("4T3")       
    elif phase_angle > 126 and phase_angle < 144:
        summary_fura_ring_conf_from_C_P.write("4T3---4E")
    elif phase_angle == 144:
        summary_fura_ring_conf_from_C_P.write("4E(C4'-endo)")   
    elif phase_angle > 144 and phase_angle < 162:
        summary_fura_ring_conf_from_C_P.write("4E---4T0")
    elif phase_angle == 162:
        summary_fura_ring_conf_from_C_P.write("4T0")         
    elif phase_angle > 162 and phase_angle < 180:
        summary_fura_ring_conf_from_C_P.write("4T0---E0")
    elif phase_angle == 180:
        summary_fura_ring_conf_from_C_P.write("E0(O4-exo)")
    elif phase_angle > 180 and phase_angle < 198:
        summary_fura_ring_conf_from_C_P.write("E0---1T0")
    elif phase_angle == 198:
        summary_fura_ring_conf_from_C_P.write("1T0")
    elif phase_angle > 198 and phase_angle < 216:
        summary_fura_ring_conf_from_C_P.write("1T0---1E")
    elif phase_angle == 216:
        summary_fura_ring_conf_from_C_P.write("1E(C1'-endo)")
    elif phase_angle > 216 and phase_angle < 234:
        summary_fura_ring_conf_from_C_P.write("1E---1T2")
    elif phase_angle == 234:
        summary_fura_ring_conf_from_C_P.write("1T2")
    elif phase_angle > 234 and phase_angle < 252:
        summary_fura_ring_conf_from_C_P.write("1T2---E2")
    elif phase_angle == 252:
        summary_fura_ring_conf_from_C_P.write("E2(C2'-exo)")
    elif phase_angle > 252 and phase_angle < 270:
        summary_fura_ring_conf_from_C_P.write("E2---3T2(North)")
    elif phase_angle == 270:
        summary_fura_ring_conf_from_C_P.write("3T2(C2'-exo_C3'-endo)(North)")  
    elif phase_angle > 270 and phase_angle < 288:
        summary_fura_ring_conf_from_C_P.write("3T2---3E(North)")
    elif phase_angle == 288:
        summary_fura_ring_conf_from_C_P.write("3E(C3'-endo)(North)")  
    elif phase_angle > 288 and phase_angle < 306:
        summary_fura_ring_conf_from_C_P.write("3E---3T4")
    elif phase_angle == 306:    
        summary_fura_ring_conf_from_C_P.write("3T4")
    elif phase_angle > 306 and phase_angle < 324: 
        summary_fura_ring_conf_from_C_P.write("3T4---E4")
    elif phase_angle == 324:    
        summary_fura_ring_conf_from_C_P.write("E4(C4'-exo)")
    elif phase_angle > 324 and phase_angle < 342:
        summary_fura_ring_conf_from_C_P.write("E4--0T4")
    elif phase_angle == 342:    
        summary_fura_ring_conf_from_C_P.write("0T4")
    elif phase_angle > 342 and phase_angle < 360:
        summary_fura_ring_conf_from_C_P.write("E2(C2'-exo)---3T2")
    elif phase_angle == 360:
        summary_fura_ring_conf_from_C_P.write("0E(O4-endo)")   


END_SCRIPT

# Run script.py with python3
python3 calc_ring_puckering.py

}






########################################################
##                                                    ##
## FUNCTION TO OBTAINS THE RING PUCKERING PARAMETERS  ##
## AND MOST PROBABLE SUGAR RING CONFORMATION.         ##
##                                                    ##
########################################################

function obtain_ring_atoms_coord {

# Load function input variables
local environment=$1
local TC=$2
local ring_conf=$3
local RU=$4

if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ]; then
# Obtain the 5-member ring atoms labels  
C1=8
C2=7
C3=5
C4=3
O4=4

# Copy cartesian coordinates of this atoms in corresponding files
sed -n ''$C1'p' cartesian_coordinates.txt > C1_fura_ring.txt
sed -n ''$C2'p' cartesian_coordinates.txt > C2_fura_ring.txt
sed -n ''$C3'p' cartesian_coordinates.txt > C3_fura_ring.txt
sed -n ''$C4'p' cartesian_coordinates.txt > C4_fura_ring.txt
sed -n ''$O4'p' cartesian_coordinates.txt > O4_fura_ring.txt

calc_puckering_param $TC

elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then
# Obtain the 6-member ring atoms labels  
C1=1
C2=2
C3=3
C4=4
C5=5
O4=6

# Copy cartesian coordinates of this atoms in corresponding files
sed -n ''$C1'p' cartesian_coordinates.txt > C1_pyra_ring.txt
sed -n ''$C2'p' cartesian_coordinates.txt > C2_pyra_ring.txt
sed -n ''$C3'p' cartesian_coordinates.txt > C3_pyra_ring.txt
sed -n ''$C4'p' cartesian_coordinates.txt > C4_pyra_ring.txt
sed -n ''$C5'p' cartesian_coordinates.txt > C5_pyra_ring.txt
sed -n ''$O4'p' cartesian_coordinates.txt > O4_pyra_ring.txt

calc_puckering_param $TC

# Copy cartesian coordinates of this atoms in corresponding files

elif [ "$TC" == "threose" ]; then
# Obtain the 5-member ring atoms labels  
C1=5
C2=4
C3=3
C4=1
O4=2

# Copy cartesian coordinates of this atoms in corresponding files
sed -n ''$C1'p' cartesian_coordinates.txt > C1_fura_ring.txt
sed -n ''$C2'p' cartesian_coordinates.txt > C2_fura_ring.txt
sed -n ''$C3'p' cartesian_coordinates.txt > C3_fura_ring.txt
sed -n ''$C4'p' cartesian_coordinates.txt > C4_fura_ring.txt
sed -n ''$O4'p' cartesian_coordinates.txt > O4_fura_ring.txt

calc_puckering_param TC

fi


}





function w_results_to_analysis_folder {

# Load function input variables
local environment=$1
local TC=$2
local ring_conf=$3
local RU=$4

# Send results to general files outside in nucleosides folder
if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

# Load Cremer & Pople parameters
ring_conformation_C_P=`cat 5-member_ring_conformation_C_P.txt`
q2_amplitude_C_P=`cat amplitude_q2_from_C_P.txt`
total_amplitude_Q_C_P=`cat total_amplitude_Q_from_C_P.txt`
phase_angle_phi2_C_P=`cat phase_angle_Phi2.txt`

# Output variables in general files inside the environment folder depending on the ring type
echo "$ring_conformation_C_P" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/final_ring_conformation_C_P_furanoses.txt
echo "$q2_amplitude_C_P" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/q2_C_P_furanoses.txt
echo "$total_amplitude_Q_C_P" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/total_amplitude_Q_C_P_furanoses.txt
echo "$phase_angle_phi2_C_P" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/phase_angle_phi2_Q_C_P_furanoses.txt
echo "$RU" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/RU_names_furanoses.txt
echo "$TC" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/TC_names_furanoses.txt
echo "$ring_conf" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/initial_ring_conf_furanoses.txt
echo "$environment" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/environment_list_furanoses.txt




elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then

# Load C & P parameters
ring_conformation_pyra_C_P=`cat 6-member_ring_conformation_C_P_from_theta.txt`
q2_amplitude_pyra_C_P=`cat amplitude_q2_from_C_P.txt`
q3_amplitude_pyra_C_P=`cat amplitude_q3_from_C_P.txt`
total_amplitude_Q_pyra_C_P=`cat total_amplitude_Q_from_C_P.txt`
phase_angle_phi2_pyra_C_P=`cat phase_angle_Phi2.txt`
phase_angle_theta_C_P=`cat phase_angle_theta.txt`

# Output variables in general files inside the environment folder depending on the ring type
echo "$ring_conformation_pyra_C_P" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/final_ring_conformation_C_P_pyranoses.txt
echo "$q2_amplitude_pyra_C_P" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/q2_C_P_pyranoses.txt
echo "$q3_amplitude_pyra_C_P" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/q3_C_P_pyranoses.txt
echo "$total_amplitude_Q_pyra_C_P" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/total_amplitude_Q_C_P_pyranoses.txt
echo "$phase_angle_phi2_pyra_C_P" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/phase_angle_phi2_Q_C_P_pyranoses.txt
echo "$phase_angle_theta_C_P" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/phase_angle_theta_Q_C_P_pyranoses.txt
echo "$RU" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/RU_names_pyranoses.txt
echo "$TC" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/TC_names_pyranoses.txt
echo "$ring_conf" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/initial_ring_conf_pyranoses.txt
echo "$environment" >> $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC/environment_list_pyranoses.txt


fi

}






##########################################################################################################################
############################################## WORK ON EACH FOLDER #######################################################
##########################################################################################################################

# Create folder for analysis and post-processing
mkdir $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering

# Open corresponding folder
cd nucleosides/final_gaussian_opt

for environment in "${environment[@]}"; do


cd $environment

for TC in "${TC_names[@]}"; do

# Create corresponding folder for furanoses
mkdir $MD_RUN/nucleosides/post_proc_summary_results_csv_graphs/sugar_puckering/$TC

cd $TC

if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

for ring_conf_fura in "${Ring_Conf_furanose[@]}"; do

cd $ring_conf_fura

for RU in "${RU_names[@]}"; do

cd $RU

# Output message of status of script
echo "working on $environment/$TC/$ring_conf_fura/$RU"

# Run function to obtain the ring atoms coordinates
obtain_ring_atoms_coord $environment $TC $ring_conf_fura $RU 

# Run function to calculate the ring puckering parameters
calc_puckering_param $environment $TC $ring_conf_fura $RU  

# Write results to general files in analysis folder
w_results_to_analysis_folder $environment $TC $ring_conf_fura $RU


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

# Run function to obtain the ring atoms coordinates
obtain_ring_atoms_coord $environment $TC $ring_conf_pyra $RU  

# Run function to calculate the ring puckering parameters
calc_puckering_param $environment $TC $ring_conf_pyra $RU  

# Write results to general files in analysis folder
w_results_to_analysis_folder $environment $TC $ring_conf_pyra $RU  


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


#####################################
###### Generate final csv files #####
#####################################

cd post_proc_summary_results_csv_graphs/sugar_puckering

for TC in "${TC_names[@]}"; do

cd $TC

# Delete the [ and ] by from some of the files
if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then
sed -i 's/[][]//g' phase_angle_phi2_Q_C_P_furanoses.txt
sed -i 's/[][]//g' q2_C_P_furanoses.txt
else
sed -i 's/[][]//g' phase_angle_phi2_Q_C_P_pyranoses.txt
fi

# Create and run python script to gen summary csv files and freq hist 
export TC

cat >gen_ring_puckering_csv_graph.py <<'END_SCRIPT'

# Import environmental variables from python libraries
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
from matplotlib import ticker, cm
from matplotlib.transforms import Bbox
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from matplotlib.ticker import PercentFormatter
from matplotlib import ticker as tick
from matplotlib.lines import Line2D
from matplotlib.gridspec import SubplotSpec
import matplotlib.lines as mlines


# Import environmental variables from bash script for angles in degrees
TC = os.environ.get("TC")
#TC = "2deoxy_ribopyranose"

# Create important functions
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
    
    
def create_subtitle(fig: plt.Figure, grid: SubplotSpec, title: str):
    "Sign sets of subplots with title"
    row = fig.add_subplot(grid)
    # the '\n' is important
    #if title == 'HAsO$_3^-$-2dRibf-RU in vacuum':
    row.set_title(f'{title}\n', fontweight='bold', fontsize=28)
    #else:
     #   row.set_title(f'{title}', fontweight='semibold', fontsize=26)
    # hide subplot
    row.set_frame_on(False)
    row.axis('off') 
    
    
def realign_polar_xticks(ax):
    for theta, label in zip(ax.get_xticks(), ax.get_xticklabels()):
        theta = theta * ax.get_theta_direction() + ax.get_theta_offset()
        theta = np.pi/2 - theta
        y, x = np.cos(theta), np.sin(theta)
        if x >= 0.1:
            label.set_horizontalalignment('left')
        if x <= -0.1:
            label.set_horizontalalignment('right')
        if y >= 0.5:
            label.set_verticalalignment('bottom')
        if y <= -0.5:
            label.set_verticalalignment('top')  




          
# Work on each TC

if TC in ["2deoxy_ribofuranose", "ribofuranose", "threose"]:

    ############## Create csv file for furanose ring data #################
    environment = np.loadtxt("environment_list_furanoses.txt", dtype='str')
    TC_fura = np.loadtxt("TC_names_furanoses.txt", dtype='str')
    initial_ring_conf_fura = np.loadtxt("initial_ring_conf_furanoses.txt", dtype='str')
    RU_fura = np.loadtxt("RU_names_furanoses.txt", dtype='str')
    final_ring_conf_fura = np.loadtxt("final_ring_conformation_C_P_furanoses.txt", dtype='str')
    q2_fura_C_P = np.loadtxt("q2_C_P_furanoses.txt", dtype='float')
    Q_fura_C_P = np.loadtxt("total_amplitude_Q_C_P_furanoses.txt", dtype='float')
    phi2_fura_C_P = np.loadtxt("phase_angle_phi2_Q_C_P_furanoses.txt", dtype='float')
    
    puckering_df = pd.DataFrame({'Environment': environment, 'TC': TC_fura, 
                                 'Initial_ring_conf': initial_ring_conf_fura, 
                                 'RU': RU_fura, 'q2(amplitude)_Cremer&Pople': q2_fura_C_P, 
                                 'Q(Total_amplitude)_Cremer&Pople': Q_fura_C_P, 
                                 'phi2(Phase_angle)_Cremer&Pople':phi2_fura_C_P, 
                                 'Final_ring_conf_Cremer&Pople': final_ring_conf_fura})

    # Convert pandas dataframe to csv file
    puckering_df.to_csv('db_'+str(TC)+'_ring_puckering_params.csv', index=False)
    
        ############### GENERATE GRAPHS ################
    
    ################ Obtain entries for vacuum and water environment from db #################
    db_vacuum = puckering_df[puckering_df['Environment'] == "vacuum"]
    db_water = puckering_df[puckering_df['Environment'] == "water"]
       
    ################# Obtain from db_env the entries for canonical and for not_canonical bases ################
    # List of RUs
    canonical_RU = ['adenine', 'guanine', 'cytosine', 'thymine', 'uracil']
    non_canonical_RU = ['TARC_cbond', 'TARC_nbond', 'BA_cbond', 'BA_nbond', 'melamine', 'CA']
    
    #vacuum
    db_vacuum_canonical = db_vacuum[db_vacuum['RU'].isin(canonical_RU)]
    db_vacuum_non_canonical = db_vacuum[db_vacuum['RU'].isin(non_canonical_RU)]
    
    # Water
    db_water_canonical = db_water[db_water['RU'].isin(canonical_RU)]
    db_water_non_canonical = db_water[db_water['RU'].isin(non_canonical_RU)]
    
    ################ Obtain from db_env the entries for beta and alpha anomers ##################
    alpha = ["2endo_alpha", "3endo_alpha"]
    beta = ["2endo_beta", "3endo_beta"]
    
    # Vacuum
    all_vacuum_alpha = db_vacuum[db_vacuum['Initial_ring_conf'].isin(alpha)]
    all_vacuum_beta = db_vacuum[db_vacuum['Initial_ring_conf'].isin(beta)]   
    can_vacuum_alpha = db_vacuum_canonical[db_vacuum_canonical['Initial_ring_conf'].isin(alpha)]
    can_vacuum_beta = db_vacuum_canonical[db_vacuum_canonical['Initial_ring_conf'].isin(beta)]
    non_can_vacuum_alpha = db_vacuum_non_canonical[db_vacuum_non_canonical['Initial_ring_conf'].isin(alpha)]
    non_can_vacuum_beta = db_vacuum_non_canonical[db_vacuum_non_canonical['Initial_ring_conf'].isin(beta)]
    
    # Water
    all_water_alpha = db_water[db_water['Initial_ring_conf'].isin(alpha)]
    all_water_beta = db_water[db_water['Initial_ring_conf'].isin(beta)]   
    can_water_alpha = db_water_canonical[db_water_canonical['Initial_ring_conf'].isin(alpha)]
    can_water_beta = db_water_canonical[db_water_canonical['Initial_ring_conf'].isin(beta)]
    non_can_water_alpha = db_water_non_canonical[db_water_non_canonical['Initial_ring_conf'].isin(alpha)]
    non_can_water_beta = db_water_non_canonical[db_water_non_canonical['Initial_ring_conf'].isin(beta)]
      
    ################### Obtain values of phi2 angle #########################
    phi_all_vacuum_alpha = all_vacuum_alpha['phi2(Phase_angle)_Cremer&Pople']
    phi_all_vacuum_beta = all_vacuum_beta['phi2(Phase_angle)_Cremer&Pople']
    phi_can_vacuum_alpha = can_vacuum_alpha ['phi2(Phase_angle)_Cremer&Pople']
    phi_can_vacuum_beta = can_vacuum_beta ['phi2(Phase_angle)_Cremer&Pople']
    phi_non_can_vacuum_alpha = non_can_vacuum_alpha ['phi2(Phase_angle)_Cremer&Pople']
    phi_non_can_vacuum_beta = non_can_vacuum_beta ['phi2(Phase_angle)_Cremer&Pople']
    
    # Water
    phi_all_water_alpha = all_water_alpha['phi2(Phase_angle)_Cremer&Pople']
    phi_all_water_beta = all_water_beta['phi2(Phase_angle)_Cremer&Pople']
    phi_can_water_alpha = can_water_alpha ['phi2(Phase_angle)_Cremer&Pople']
    phi_can_water_beta = can_water_beta ['phi2(Phase_angle)_Cremer&Pople']
    phi_non_can_water_alpha = non_can_water_alpha ['phi2(Phase_angle)_Cremer&Pople']
    phi_non_can_water_beta = non_can_water_beta ['phi2(Phase_angle)_Cremer&Pople'] 
    
    ############### Obtain the frequency tables and counts ##################
    # Vacuum
    (freq_df_all_vacuum_alpha, count_phi_all_vacuum_alpha) = frequencies_for_hist(phi_all_vacuum_alpha, 
                                                                                  0, 360, 18)
    (freq_df_all_vacuum_beta, count_phi_all_vacuum_beta) = frequencies_for_hist(phi_all_vacuum_beta, 
                                                                                  0, 360, 18)  
                                                                                                                                                                                                                                                                                                                                           
    (freq_df_can_vacuum_alpha, count_phi_can_vacuum_alpha) = frequencies_for_hist(phi_can_vacuum_alpha, 
                                                                                  0, 360, 18)
    (freq_df_can_vacuum_beta, count_phi_can_vacuum_beta) = frequencies_for_hist(phi_can_vacuum_beta, 
                                                                                  0, 360, 18)                                                                             
    (freq_df_non_can_vacuum_alpha, count_phi_non_canonical_vacuum_alpha) = frequencies_for_hist(phi_non_can_vacuum_alpha,
                                                                                                0, 360, 18)
    (freq_df_non_can_vacuum_beta, count_phi_non_canonical_vacuum_beta) = frequencies_for_hist(phi_non_can_vacuum_beta,
                                                                                                0, 360, 18) 
                                                                                                                                                                                             
    # Water
    (freq_df_all_water_alpha, count_phi_all_water_alpha) = frequencies_for_hist(phi_all_water_alpha, 
                                                                                0, 360, 18)
    (freq_df_all_water_beta, count_phi_all_water_beta) = frequencies_for_hist(phi_all_water_beta, 
                                                                              0, 360, 18)  
                                                                                                                                                                                                                                                                                                                                           
    (freq_df_can_water_alpha, count_phi_can_water_alpha) = frequencies_for_hist(phi_can_water_alpha, 
                                                                                0, 360, 18)
    (freq_df_can_water_beta, count_phi_can_water_beta) = frequencies_for_hist(phi_can_water_beta, 
                                                                              0, 360, 18)                                                                             
    (freq_df_non_can_water_alpha, count_phi_non_canonical_water_alpha) = frequencies_for_hist(phi_non_can_water_alpha,
                                                                                              0, 360, 18)
    (freq_df_non_can_water_beta, count_phi_non_canonical_water_beta) = frequencies_for_hist(phi_non_can_water_beta,
                                                                                            0, 360, 18) 
    
    ################## Obtain the circ mean, std and total datapoints #####################
    # Vacuum   
    phi_can_vacuum_rad = np.radians(db_vacuum_canonical['phi2(Phase_angle)_Cremer&Pople'])
    (circ_mean_vacuum_canonical, circ_std_vacuum_canonical) = (np.degrees(circmean(phi_can_vacuum_rad)),
                                                               np.degrees(circstd(phi_can_vacuum_rad)))
                                                                    
    phi_non_can_vacuum_rad = np.radians(db_vacuum_non_canonical['phi2(Phase_angle)_Cremer&Pople'])    
    (circ_mean_vacuum_non_canonical, circ_std_vacuum_non_canonical) = (np.degrees(circmean(phi_non_can_vacuum_rad)),
                                                                       np.degrees(circstd(phi_non_can_vacuum_rad)))

    phi_all_vacuum_rad = np.radians(db_vacuum['phi2(Phase_angle)_Cremer&Pople'])    
    (circ_mean_all_vacuum, circ_std_all_vacuum) = (np.degrees(circmean(phi_all_vacuum_rad)),
                                                                       np.degrees(circstd(phi_all_vacuum_rad)))

    # Water
    phi_can_water_rad = np.radians(db_water_canonical['phi2(Phase_angle)_Cremer&Pople'])
    (circ_mean_water_canonical, circ_std_water_canonical) = (np.degrees(circmean(phi_can_water_rad)),
                                                             np.degrees(circstd(phi_can_water_rad)))
                                                                      
    phi_non_can_water_rad = np.radians(db_water_non_canonical['phi2(Phase_angle)_Cremer&Pople'])    
    (circ_mean_water_non_canonical, circ_std_water_non_canonical) = (np.degrees(circmean(phi_non_can_water_rad)),
                                                                     np.degrees(circstd(phi_non_can_water_rad)))

    phi_all_water_rad = np.radians(db_water['phi2(Phase_angle)_Cremer&Pople'])    
    (circ_mean_all_water, circ_std_all_water) = (np.degrees(circmean(phi_all_water_rad)),
                                                                     np.degrees(circstd(phi_all_water_rad)))
                                                                     
    info_can_vacuum = [count_phi_can_vacuum_alpha, count_phi_can_vacuum_beta, circ_mean_vacuum_canonical, circ_std_vacuum_canonical]
    info_non_can_vacuum = [count_phi_non_canonical_vacuum_alpha, count_phi_non_canonical_vacuum_beta, circ_mean_vacuum_non_canonical, circ_std_vacuum_non_canonical]
    info_all_vacuum = [count_phi_all_vacuum_alpha, count_phi_all_vacuum_beta, circ_mean_all_vacuum, circ_std_all_vacuum]
    
    info_can_water = [count_phi_can_water_alpha, count_phi_can_water_beta, circ_mean_water_canonical, circ_std_water_canonical]
    info_non_can_water = [count_phi_non_canonical_water_alpha, count_phi_non_canonical_water_beta, circ_mean_water_non_canonical, circ_std_water_non_canonical]
    info_all_water = [count_phi_all_water_alpha, count_phi_all_water_beta, circ_mean_all_water, circ_std_all_water]

    ############# Generate rose diagram ###############
    
    # Generate the graphs
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
        theta = np.radians(np.arange(0, 360, 18)) 
        mean = round(mean, 1)
        std = round(std, 1)
        
        # Obtain the range for freqs
        width = np.radians(18)
    
        # Create bar plot
        bars = ax.bar(theta, count_alpha, width = width, edgecolor = 'black', align = 'edge', facecolor='gray', 
                      linewidth = 2.0, label = 'Alpha')
            
        bars2 = ax.bar(theta, count_beta, width = width, edgecolor = 'black', align = 'edge', facecolor='g', 
                       bottom = count_alpha, linewidth = 2.0, label = 'Beta')
                       
        # Modify grid 
        #ax.grid(zorder=0)
        ax.grid(color='grey', linestyle='--')
    
        labels=['0°', "18°", '36°', "54°", '72°', "90°", '108°', "126°", '144°', "162°", 
         "180°", "198°", '216°', "234°", '252°', "270°", '288°', "306°", '324°', "342°"]
        ax.set_xticks(theta)
        #ax.set_xticklabels(labels, fontweight='bold')
        ax.set_xticklabels(labels)
        
        # Create red arrow that points the circular mean 
        x_tail = 0.0
        y_tail = 0.0
        x_head = float(mean)/180.*np.pi
        array = count_alpha + count_beta
        y_head = max(array)
        dx = x_head - x_tail
        dy = y_head - y_tail
        arrow = mpatches.FancyArrowPatch((x_tail, y_tail), (dx, dy),
                                          mutation_scale=30, facecolor = 'red', edgecolor='black')
        ax.add_patch(arrow)
        
        # Add the title to rose diagram
        if i == 0:
            ax.set_title('Canonical RUs' + '\n\n' + '(' + r'$\bar\phi_2$' + '=' + str(mean) + "°" +  "±" + str(std) + "°" +
                         ")" + '\n', fontweight="bold", fontsize=28)
        elif i == 1:
            ax.set_title('\nNon-canonical RUs' + '\n\n' + '(' + r'$\bar\phi_2$' + '=' + str(mean) + "°" + "±" + str(std) + "°" +
                         ")" + '\n', fontweight="bold", fontsize=28)
        elif i == 2:
            ax.set_title('\nAll RUs' + '\n\n' + '(' + r'$\bar\phi_2$' + '=' + str(mean) + "°" + "±" + str(std) + "°" + ")" + 
                         '\n', fontweight="bold", fontsize=28)
        elif i == 3:
            ax.set_title('\n(' + r'$\bar\phi_2$' + '=' + str(mean) + "°" + "±" + str(std) + "°" + ")" +
                         '\n', fontweight="bold", fontsize=28)
        elif i == 4:
            ax.set_title('\n(' + r'$\bar\phi_2$' + '=' + str(mean) + "°" + "±" + str(std) + "°" + ")" 
                         + '\n', fontweight="bold", fontsize=28)
        elif i == 5:
            ax.set_title('\n(' + r'$\bar\phi_2$' + '=' + str(mean) + "°" + "±" + str(std) + "°" + ")" +
                         '\n', fontweight="bold", fontsize=28)
            #ax.legend(fontsize=24, bbox_to_anchor=(1.55, 0.00), loc="lower right")
            ax.legend(['α', 'β'], fontsize=26, bbox_to_anchor=(1.5, 0.00), loc="lower right")
            leg = ax.get_legend()
            leg.legendHandles[0].set_color('grey')
            leg.legendHandles[1].set_color('green') 
            
        # little space between the bar and the label
        labelPadding = 4
        angles = np.radians(np.arange(9, 360, 18)) 
        
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
                    fontsize=22,
                    fontweight='bold')
        
        #ax.set_theta_direction(-1)
                    
        # Set the bar numbers (First way)
        for bar, angle, height, label, ad_height in zip(bars2, angles, count_beta, count_beta, count_alpha):        
            # Labels are rotated. Rotation must be specified in degrees :(
            rotation = np.rad2deg(angle)        
            # Flip some labels upside down
            alignment = ""
            if angle >= np.pi/2 and angle < 3*np.pi/2:
                alignment = "right"
                rotation = rotation - 180
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
                    fontsize=22,
                    fontweight='bold')
                    
        # Create space between circle and angles vals           
        #realign_polar_xticks (ax)
        
        # Configure the axis ticks
        ax.tick_params(axis='x', colors='black', labelsize='26')
        ax.tick_params(axis='y', direction='out', colors='blue', labelsize='22')
        #plt.xticks(weight = 'bold')
        labelss = ax.get_yticklabels()
        for label in labelss:
            label.set_fontweight('bold') 
            
        # Change position of y labels
        ax.set_rlabel_position(90)    
            
        # Top of the the bars you can do the following:
        for k, spine in ax.spines.items():  #ax.spines is a dictionary
            spine.set_zorder(15) 
        # Change position of y labels
        #if i == 5:
        #ax.set_rlabel_position(270)
        #else:
        #ax.set_rlabel_position(90)
                                           
    plt.tight_layout()    
    plt.savefig('rose_diagram_'+str(TC)+'.pdf')                                                                                 

else: # Case of pyranose rings

    ############## Create csv file for pyranose ring data #################
    environment = np.loadtxt("environment_list_pyranoses.txt", dtype='str')
    TC_pyra = np.loadtxt("TC_names_pyranoses.txt", dtype='str')
    initial_ring_conf_pyra = np.loadtxt("initial_ring_conf_pyranoses.txt", dtype='str')
    RU_pyra = np.loadtxt("RU_names_pyranoses.txt", dtype='str')
    final_ring_conf_pyra = np.loadtxt("final_ring_conformation_C_P_pyranoses.txt", dtype='str')    
    q2_pyra_C_P = np.loadtxt("q2_C_P_pyranoses.txt", dtype='float')
    q3_pyra_C_P = np.loadtxt("q3_C_P_pyranoses.txt", dtype='float')
    Q_pyra_C_P = np.loadtxt("total_amplitude_Q_C_P_pyranoses.txt", dtype='float')
    phi2_pyra_C_P = np.loadtxt("phase_angle_phi2_Q_C_P_pyranoses.txt", dtype='str')
    theta_pyra_C_P = np.loadtxt("phase_angle_theta_Q_C_P_pyranoses.txt", dtype='float')

    puckering_df = pd.DataFrame({'Environment':environment, 'TC': TC_pyra, 
                                 'Initial_ring_conf': initial_ring_conf_pyra, 
                                 'RU': RU_pyra, 
                                 'q2(amplitude)_Cremer&Pople': q2_pyra_C_P, 
                                 'q3(amplitude)_Cremer&Pople': q3_pyra_C_P, 
                                 'Q(Total_amplitude)_Cremer&Pople': Q_pyra_C_P, 
                                 'phi2(phase_angle)_Cremer&Pople':phi2_pyra_C_P, 
                                 'theta(phase_angle)_Cremer&Pople':theta_pyra_C_P, 
                                 'Final_ring_conf_Cremer&Pople': final_ring_conf_pyra})
    # Convert pandas dataframe to csv file
    puckering_df.to_csv('db_'+str(TC)+'_ring_puckering_params.csv', index=False)
    
    ################ Obtain entries for vacuum and water environment from db #################
    db_vacuum = puckering_df[puckering_df['Environment'] == "vacuum"]
    db_water = puckering_df[puckering_df['Environment'] == "water"]
    
    ################# Obtain from db_env the entries for canonical and for not_canonical bases ################
    # List of RUs
    canonical_RU = ['adenine', 'guanine', 'cytosine', 'thymine', 'uracil']
    non_canonical_RU = ['TARC_cbond', 'TARC_nbond', 'BA_cbond', 'BA_nbond', 'melamine', 'CA']
    
    #vacuum
    db_vacuum_canonical = db_vacuum[db_vacuum['RU'].isin(canonical_RU)]
    db_vacuum_non_canonical = db_vacuum[db_vacuum['RU'].isin(non_canonical_RU)]
    
    # Water
    db_water_canonical = db_water[db_water['RU'].isin(canonical_RU)]
    db_water_non_canonical = db_water[db_water['RU'].isin(non_canonical_RU)]
    
    ################# Obtain cells for alpha and beta anomers ################
    ################ Obtain from db_env the entries for beta and alpha anomers ##################
    alpha = ["D1C4_alpha", "D4C1_alpha"]
    beta = ["D1C4_beta", "D4C1_beta"]
    
    # Vacuum
    db_vacuum_canonical_alpha = db_vacuum_canonical[db_vacuum_canonical['Initial_ring_conf'].isin(alpha)]
    db_vacuum_canonical_beta = db_vacuum_canonical[db_vacuum_canonical['Initial_ring_conf'].isin(beta)] 
    db_vacuum_non_canonical_alpha = db_vacuum_non_canonical[db_vacuum_non_canonical['Initial_ring_conf'].isin(alpha)]
    db_vacuum_non_canonical_beta = db_vacuum_non_canonical[db_vacuum_non_canonical['Initial_ring_conf'].isin(beta)]
    db_vacuum_alpha = db_vacuum[db_vacuum['Initial_ring_conf'].isin(alpha)]
    db_vacuum_beta = db_vacuum[db_vacuum['Initial_ring_conf'].isin(beta)]

    # Water
    db_water_canonical_alpha = db_water_canonical[db_water_canonical['Initial_ring_conf'].isin(alpha)]
    db_water_canonical_beta = db_water_canonical[db_water_canonical['Initial_ring_conf'].isin(beta)] 
    db_water_non_canonical_alpha = db_water_non_canonical[db_water_non_canonical['Initial_ring_conf'].isin(alpha)]
    db_water_non_canonical_beta = db_water_non_canonical[db_water_non_canonical['Initial_ring_conf'].isin(beta)]
    db_water_alpha = db_water[db_water['Initial_ring_conf'].isin(alpha)]
    db_water_beta = db_water[db_water['Initial_ring_conf'].isin(beta)]
    
    ############### Generate the heatmaps ##########
    ################### Obtain values of phi2 and theta angles #########################
    phi2_vacuum_can = db_vacuum_canonical['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_vacuum_can = db_vacuum_canonical['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_vacuum_non_can = db_vacuum_non_canonical['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_vacuum_non_can = db_vacuum_non_canonical['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_vacuum_all = db_vacuum['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_vacuum_all = db_vacuum['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    
    phi2_water_can = db_water_canonical['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_water_can = db_water_canonical['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_water_non_can = db_water_non_canonical['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_water_non_can = db_water_non_canonical['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_water_all = db_water['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_water_all = db_water['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    
    ################### Obtain values of phi2 and theta angles in alpha and beta #########################
    phi2_vacuum_can_alpha = db_vacuum_canonical_alpha['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_vacuum_can_beta = db_vacuum_canonical_beta['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_vacuum_can_alpha = db_vacuum_canonical_alpha['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_vacuum_can_beta = db_vacuum_canonical_beta['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_vacuum_non_can_alpha = db_vacuum_non_canonical_alpha['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_vacuum_non_can_beta = db_vacuum_non_canonical_beta['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_vacuum_non_can_alpha = db_vacuum_non_canonical_alpha['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_vacuum_non_can_beta = db_vacuum_non_canonical_beta['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_vacuum_all_alpha = db_vacuum_alpha['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_vacuum_all_beta = db_vacuum_beta['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_vacuum_all_alpha = db_vacuum_alpha['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_vacuum_all_beta = db_vacuum_beta['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    
    phi2_water_can_alpha = db_water_canonical_alpha['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_water_can_beta = db_water_canonical_beta['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_water_can_alpha = db_water_canonical_alpha['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_water_can_beta = db_water_canonical_beta['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_water_non_can_alpha = db_water_non_canonical_alpha['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_water_non_can_beta = db_water_non_canonical_beta['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_water_non_can_alpha = db_water_non_canonical_alpha['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_water_non_can_beta = db_water_non_canonical_beta['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_water_all_alpha = db_water_alpha['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    phi2_water_all_beta = db_water_beta['phi2(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_water_all_alpha = db_water_alpha['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    theta_water_all_beta = db_water_beta['theta(phase_angle)_Cremer&Pople'].to_numpy(dtype = 'float')
    
    class RemainderFixed(axes_size.Scaled):
        def __init__(self, xsizes, ysizes, divider):
            self.xsizes =xsizes
            self.ysizes =ysizes
            self.div = divider

        def get_size(self, renderer):
            xrel, xabs = axes_size.AddList(self.xsizes).get_size(renderer)
            yrel, yabs = axes_size.AddList(self.ysizes).get_size(renderer)
            bb = Bbox.from_bounds(*self.div.get_position()).transformed(self.div._fig.transFigure)
            w = bb.width/self.div._fig.dpi - xabs
            h = bb.height/self.div._fig.dpi - yabs
            return 0, min([w,h])


    def make_square_axes_with_colorbar(ax, size=0.1, pad=0.1):
        """ Make an axes square, add a colorbar axes next to it, 
        Parameters: size: Size of colorbar axes in inches
                    pad : Padding between axes and cbar in inches
        Returns: colorbar axes
        """
        divider = make_axes_locatable(ax)
        margin_size = axes_size.Fixed(size)
        pad_size = axes_size.Fixed(pad)
        xsizes = [pad_size, margin_size]
        yhax = divider.append_axes("right", size=margin_size, pad=pad_size)
        divider.set_horizontal([RemainderFixed(xsizes, [], divider)] + xsizes)
        divider.set_vertical([RemainderFixed(xsizes, [], divider)])
        return yhax                                                                            
                                                                                       
    rows = 2
    cols = 3                                                                           
    fig, axn = plt.subplots(rows, cols, sharey=False, sharex=False, figsize=(36,22))
    grid = plt.GridSpec(rows, cols)
    #if TC == "2deoxy_ribopyranose":
        #create_subtitle(fig, grid[0, ::], 'HPO$_3^-$-2dRib-RU (a+b) in vacuum\n\n')
     #   create_subtitle(fig, grid[1, ::], 'HPO$_3^-$-2dRib-RU (a+b) in water\n')
   # elif TC == "ribopyranose":
    #    create_subtitle(fig, grid[0, ::], 'HPO$_3^-$-Rib-RU (a+b) in vacuum\n\n')
     #   create_subtitle(fig, grid[1, ::], 'HPO$_3^-$-Rib-RU (a+b) in water\n')
      
    for i, (ax, phi2, theta, environment, phi2_alpha, phi2_beta, theta_alpha, theta_beta) in enumerate([(axn.flat[0], 
                                                                                                         phi2_vacuum_can, 
                                                                                                         theta_vacuum_can, 
                                                                                                         "vacuum", 
                                                                                                         phi2_vacuum_can_alpha, 
                                                                                                         phi2_vacuum_can_beta, 
                                                                                                         theta_vacuum_can_alpha, 
                                                                                                         theta_vacuum_can_beta),
                                                                                                        (axn.flat[1], 
                                                                                                         phi2_vacuum_non_can, 
                                                                                                         theta_vacuum_non_can, 
                                                                                                         "vacuum", 
                                                                                                         phi2_vacuum_non_can_alpha, 
                                                                                                         phi2_vacuum_non_can_beta, 
                                                                                                         theta_vacuum_non_can_alpha, 
                                                                                                         theta_vacuum_non_can_beta),
                                                                                                        (axn.flat[2], 
                                                                                                         phi2_vacuum_all, 
                                                                                                         theta_vacuum_all, 
                                                                                                         "vacuum", 
                                                                                                         phi2_vacuum_all_alpha, 
                                                                                                         phi2_vacuum_all_beta, 
                                                                                                         theta_vacuum_all_alpha, 
                                                                                                         theta_vacuum_all_beta),
                                                                                                        (axn.flat[3], 
                                                                                                         phi2_water_can, 
                                                                                                         theta_water_can, 
                                                                                                         "water", 
                                                                                                         phi2_water_can_alpha, 
                                                                                                         phi2_water_can_beta, 
                                                                                                         theta_water_can_alpha, 
                                                                                                         theta_water_can_beta),
                                                                                                        (axn.flat[4], 
                                                                                                         phi2_water_non_can, 
                                                                                                         theta_water_non_can, 
                                                                                                         "water", 
                                                                                                         phi2_water_non_can_alpha, 
                                                                                                         phi2_water_non_can_beta, 
                                                                                                         theta_water_non_can_alpha, 
                                                                                                         theta_water_non_can_beta),
                                                                                                        (axn.flat[5], 
                                                                                                         phi2_water_all, 
                                                                                                         theta_water_all, 
                                                                                                         "water", 
                                                                                                         phi2_water_all_alpha, 
                                                                                                         phi2_water_all_beta, 
                                                                                                         theta_water_all_alpha, 
                                                                                                         theta_water_all_beta)]):   
    
        # Create vectors with predeterm phi2 and theta for confs
        theta_pred = [54.7, 50.8, 54.7, 50.8, 54.7, 50.8, 54.7, 50.8, 54.7, 50.8, 54.7, 50.8, 54.7] 
        theta_pred2 = [90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90]
        theta_pred3 = [125.3, 129.2, 125.3, 129.2, 125.3, 129.2, 125.3, 129.2, 125.3, 129.2, 125.3, 129.2, 125.3]
        phi2_pred = [-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180]
        #phi2_pred2 = [-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180] 
        #phi2_pred3 = [-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180]
        #phi2_pred = [0, 30, 60, 90, 120, 150, 180, 210, 230, 260, 290, 310, 330, 360]
        
        conf_names = ['E\u2083', '\u2074H\u2083', 
                      '\u2074E', '\u2074H\u2085',
                      'E\u2085', '\u2070H\u2085',
                      '\u2070E', '\u2070H\u2081',
                      'E\u2081', '\u00b2H\u2081',
                      '\u00b2E', '\u00B2H\u2083',
                      'E\u2083'] # End of first equatorial
                      
        conf_names2 = [r'B$_{3,0}$', '\u00B9S\u2083', r'$^{1,4}$B', 
                       '\u00B9S\u2085', r'B$_{2,5}$',
                       '\u2070S\u2082', r'$^{3,0}$B', 
                       '\u00B3S\u2081', r'B$_{1,4}$',
                       '\u2075S\u2081', r'$^{2,5}$B',
                       '\u00B2S\u2080', r'B$_{3,0}$'] # End of boats conf
        
        conf_names3 = ['E\u2080', '\u00B9H\u2080',
                      '\u00b9E', '\u00b9H\u2082',
                      'E\u2082', '\u00b3H\u2082',
                      '\u00b3E', '\u00B3H\u2084',
                      'E\u2084', '\u2075H\u2084',
                      '\u2075E', '\u2075H\u2080',
                      'E\u2080']
                    
        # Add scatter for predetermined conformations
        ax.scatter (phi2_pred, theta_pred2, color='black', linewidth=1, s=80)           
        for i, type in enumerate(conf_names):
            phi = phi2_pred[i]
            the = theta_pred[i]
            ax.scatter(phi, the, color='black', linewidth=1, s=120)
            if (i % 2) == 0:
                ax.text(phi-10.5, the+3.5, type, fontsize=30, fontweight='bold', color='black')
            else:
                ax.text(phi-10.0, the-10.0, type, fontsize=30, fontweight='bold', color='black')
                
        for i, type in enumerate(conf_names2):
            phi = phi2_pred[i]
            the = theta_pred2[i]
            ax.scatter(phi, the, color='black', linewidth=1, s=80)
            if (i % 2) == 0:
                ax.text(phi-18.0, the-12.0, type, fontsize=30, fontweight='bold', color='black')
            else:
                ax.text(phi-6.5, the+3.5, type, fontsize=30, fontweight='bold', color='black')  
        
        for i, type in enumerate(conf_names3):
            phi = phi2_pred[i]
            the = theta_pred3[i]
            ax.scatter(phi, the, color='black', linewidth=1, s=80)
            if (i % 2) == 0:
                ax.text(phi-10.0, the-10.0, type, fontsize=30, fontweight='bold', color='black')
            else:
                ax.text(phi-6.5, the+3.5, type, fontsize=30, fontweight='bold', color='black')  
                        
        # Get total number of data
        N = phi2.shape[0]
    
        # Peform the kernel density estimate
        #phi2 = np.array(phi2)
        #theta = np.array(theta)
        #xx, yy = np.mgrid[min(phi2):max(phi2):(N * 1j), min(theta):max(theta):(N * 1j)]
        xx, yy = np.mgrid[-180:180:((N+2) * 1j), 0:180:(N+2 * 1j)]
        #xx, yy = np.mgrid[phi2, theta]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([phi2, theta])
        kernel = st.gaussian_kde(values, bw_method=0.4)
        f = np.reshape(kernel(positions).T, xx.shape)
    
        ## Or kernel density estimate plot instead of the contourf plot
        #image = ax.imshow(np.rot90(f), extent=[min(phi2)+5, max(phi2)+5, min(theta)+5, max(theta)+5])    
        # Contourf plot
        levels = np.linspace(f.min(), f.max(), 20)
        #cfset = ax.contourf(xx, yy, f, levels = levels, cmap = plt.cm.jet, zorder=0)
        #cfset = ax.contourf(xx, yy, f, levels = levels, cmap = 'viridis', zorder=0)
        cfset = ax.contourf(xx, yy, f, levels = levels, zorder=0, cmap = 'Blues')
        #cfset = ax.contour(xx, yy, f, levels = levels, zorder=0)
        
        # Add scatter plot
        #scatter = ax.scatter(phi2, theta, color="red", linewidth=1, edgecolor="black", s=160, marker='x')
        # Add scatter plot
        scatter = ax.scatter(phi2_alpha, theta_alpha, color="gray", linewidth=0.5, edgecolor="black", s=420, 
                             marker = 's', label = "alpha")
        # Add scatter plot
        scatter = ax.scatter(phi2_beta, theta_beta, color="green", linewidth=0.5, edgecolor="black", s=350, 
                             marker = '^', label = "beta")
                                                        
        # Label axes
        #ax.clabel(image, inline=1, fontsize=10)
        if ax == axn.flat[4] or ax == axn.flat[5]:        
            ax.set_xlabel("\n" + r'$\varphi_2$' + "(°)", fontsize = 40, fontweight='bold')
        elif ax == axn.flat[0]:
            ax.set_ylabel(r'$\theta$' + "(°)", fontsize = 40, fontweight='bold')   
        elif ax == axn.flat[3]:
            ax.set_xlabel("\n" + r'$\varphi_2$' + "(°)", fontsize = 40, fontweight='bold')
            ax.set_ylabel(r'$\theta$' + "(°)", fontsize = 40, fontweight='bold') 
            
        # Set axis ticks labels
        ax.xaxis.set_ticks(np.arange(-180, 220, 60))
        ax.yaxis.set_ticks(np.arange(0, 225, 45))
        ax.tick_params(axis='both', which='major', labelsize=35)
        
        # Adding top conformation 
        #ax.text(-20.0, 180.0+3.0, '\n\u2074C\u2081', fontsize=34, fontweight='bold', color='black')
        #ax.text(-20.0, -23.0, '\n\u00b9C\u2084', fontsize=34, fontweight='bold', color='black')
        ax.text(-20.0, 180.0+3.0, '\u00b9C\u2084', fontsize=36, color='black', fontweight='bold')
        ax.text(-20.0, -23.0, '\n\u2074C\u2081', fontsize=36, color='black', fontweight='bold')
        
        # Add grid lines
        #ax.xaxis.grid(True)
        #ax.yaxis.grid(True)
        ax.grid(which='major', linewidth=0.8, color='black', linestyle='--', alpha = 0.8)
        
        # Add titles
        if ax == axn.flat[0]:
            ax.set_title('Canonical RUs\n', fontweight="bold", fontsize=42)
        elif ax == axn.flat[1]:
            ax.set_title('Non-canonical RUs\n', fontweight="bold", fontsize=42)
        elif ax == axn.flat[2]:
            ax.set_title('All RUs\n', fontweight="bold", fontsize=42)
        #elif ax == axn.flat[3]:
         #   ax.set_title('\nCanonical, water\n', fontweight="bold", fontsize=40)
        #elif ax == axn.flat[4]:
         #   ax.set_title('\nNon canonical, water\n', fontweight="bold", fontsize=40)
        #elif ax == axn.flat[5]:
         #   ax.set_title('\nAll RUs, water\n', fontweight="bold", fontsize=40)
            
        # Add color bar  
        def fmt(x, pos):
            a, b = '{:.1e}'.format(x).split('e')
            b = int(b)
            return r'${} \times 10^{{{}}}$'.format(a, b)
        cbar = plt.colorbar(cfset, ax=ax, format=ticker.FuncFormatter(fmt))
        cbar.ax.tick_params(labelsize=24)
        #cbar.ax.set_xlabel('KERNEL')
        
        # Set legend
        if ax == axn.flat[5]:
            grey_square = mlines.Line2D([], [], color='grey', marker='s', linestyle='None',
                          markersize=20, label='α')
            green_triangle = mlines.Line2D([], [], color='green', marker='^', linestyle='None',
                             markersize=20, label='β')
            ax.legend(handles=[grey_square, green_triangle], fontsize=40, bbox_to_anchor=(1.2, -0.25), loc="lower center", ncol=2)
         
        # Invert axis
        #ax.set_ylim(ax.get_ylim()[::-1])  
        
        # Adding top conformation to y axis inverted 
        #ax.text(-20.0, 204.0, '\u00b9C\u2084', fontsize=34, fontweight='bold', color='black')
        #ax.text(-20.0, -5.0, '\n\u2074C\u2081', fontsize=34, fontweight='bold', color='black')                                                          
    
    plt.tight_layout(h_pad=1)
    plt.savefig('kernel_pdf_countour_'+str(TC)+'4.pdf')     
    
    fig, axn = plt.subplots(2, 3, sharey=False, sharex=False, figsize=(32,20))
    for i, (ax, phi2_alpha, phi2_beta, theta_alpha, theta_beta, environment) in enumerate([(axn.flat[0], 
                                                                                       phi2_vacuum_can_alpha,
                                                                                       phi2_vacuum_can_beta,
                                                                                       theta_vacuum_can_alpha,
                                                                                       theta_vacuum_can_beta,
                                                                                       "vacuum"),
                                                                                      (axn.flat[1], 
                                                                                       phi2_vacuum_non_can_alpha, 
                                                                                       phi2_vacuum_non_can_beta, 
                                                                                       theta_vacuum_non_can_alpha,
                                                                                       theta_vacuum_non_can_beta,
                                                                                       "vacuum"),
                                                                                      (axn.flat[2], 
                                                                                       phi2_vacuum_all_alpha,
                                                                                       phi2_vacuum_all_beta,
                                                                                       theta_vacuum_all_alpha,
                                                                                       theta_vacuum_all_beta,
                                                                                       "vacuum"),
                                                                                      (axn.flat[3], 
                                                                                       phi2_water_can_alpha,
                                                                                       phi2_water_can_beta,
                                                                                       theta_water_can_alpha,
                                                                                       theta_water_can_beta,
                                                                                       "water"),
                                                                                      (axn.flat[4], 
                                                                                       phi2_water_non_can_alpha, 
                                                                                       phi2_water_non_can_beta, 
                                                                                       theta_water_non_can_alpha,
                                                                                       theta_water_non_can_beta,
                                                                                       "vacuum"),
                                                                                      (axn.flat[5], 
                                                                                       phi2_water_all_alpha,
                                                                                       phi2_water_all_beta,
                                                                                       theta_water_all_alpha,
                                                                                       theta_water_all_beta,
                                                                                       "vacuum")
                                                                                     ]):
        # Create vectors with predeterm phi2 and theta for confs
        theta_pred = [54.7, 50.8, 54.7, 50.8, 54.7, 50.8, 54.7, 50.8, 54.7, 50.8, 54.7, 50.8, 54.7] 
        theta_pred2 = [90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90]
        theta_pred3 = [125.3, 129.2, 125.3, 129.2, 125.3, 129.2, 125.3, 129.2, 125.3, 129.2, 125.3, 129.2, 125.3]
        phi2_pred = [-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180]
        phi2_pred2 = [-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180] 
        phi2_pred3 = [-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180]
        
        conf_names = ['E\u2083', '\u2074H\u2083', 
                      '\u2074E', '\u2074H\u2085',
                      'E\u2085', '\u2070H\u2085',
                      '\u2070E', '\u2070H\u2081',
                      'E\u2081', '\u00b2H\u2081',
                      '\u00b2E', '\u00B2H\u2083',
                      'E\u2083'] # End of first equatorial
                      
        conf_names2 = [r'B$_{3,0}$', '\u00B9S\u2083', r'$^{1,4}$B', 
                       '\u00B9S\u2085', r'B$_{2,5}$',
                       '\u2070S\u2082', r'$^{3,0}$B', 
                       '\u00B3S\u2081', r'B$_{1,4}$',
                       '\u2075S\u2081', r'$^{2,5}$B',
                       '\u00B2S\u2080', r'B$_{3,0}$'] # End of boats conf
        
        conf_names3 = ['E\u2080', '\u00B9H\u2080',
                      '\u00b9E', '\u00b9H\u2082',
                      'E\u2082', '\u00b3H\u2082',
                      '\u00b3E', '\u00B3H\u2084',
                      'E\u2084', '\u2075H\u2084',
                      '\u2075E', '\u2075H\u2080',
                      'E\u2080']
                    
        # Add scatter for predetermined conformations
        ax.scatter (phi2_pred2, theta_pred2, color='black', linewidth=1, s=80)           
        for i, type in enumerate(conf_names):
            phi = phi2_pred[i]
            the = theta_pred[i]
            ax.scatter(phi, the, color='black', linewidth=1, s=120)
            if (i % 2) == 0:
                ax.text(phi-10.5, the+3.5, type, fontsize=30, fontweight='bold', color='black')
            else:
                ax.text(phi-10.0, the-10.0, type, fontsize=30, fontweight='bold', color='black')
                
        for i, type in enumerate(conf_names2):
            phi = phi2_pred2[i]
            the = theta_pred2[i]
            ax.scatter(phi, the, color='black', linewidth=1, s=80)
            if (i % 2) == 0:
                ax.text(phi-18.0, the-12.0, type, fontsize=30, fontweight='bold', color='black')
            else:
                ax.text(phi-6.5, the+3.5, type, fontsize=30, fontweight='bold', color='black')  
        
        for i, type in enumerate(conf_names3):
            phi = phi2_pred3[i]
            the = theta_pred3[i]
            ax.scatter(phi, the, color='black', linewidth=1, s=80)
            if (i % 2) == 0:
                ax.text(phi-10.0, the-10.0, type, fontsize=30, fontweight='bold', color='black')
            else:
                ax.text(phi-6.5, the+3.5, type, fontsize=30, fontweight='bold', color='black')  
    
        # Add scatter plot
        scatter = ax.scatter(phi2_alpha, theta_alpha, color="gray", linewidth=0.5, edgecolor="black", s=400, 
                             marker = 's', label = "alpha")
        # Add scatter plot
        scatter = ax.scatter(phi2_beta, theta_beta, color="green", linewidth=0.5, edgecolor="black", s=320, 
                             marker = '^', label = "beta")
    
        # Label axes
        #ax.clabel(image, inline=1, fontsize=10)
        if ax == axn.flat[4] or ax == axn.flat[5]:        
            ax.set_xlabel("\n" + r'$\phi$' + "(°)", fontsize = 38, fontweight='bold')
        elif ax == axn.flat[0]:
            ax.set_ylabel(r'$\theta$' + "(°)", fontsize = 38, fontweight='bold')   
        elif ax == axn.flat[3]:
            ax.set_xlabel("\n" + r'$\phi$' + "(°)", fontsize = 38, fontweight='bold')
            ax.set_ylabel(r'$\theta$' + "(°)", fontsize = 38, fontweight='bold')  
    
        # Set axis ticks labels
        ax.xaxis.set_ticks(np.arange(-180, 240, 60))
        ax.yaxis.set_ticks(np.arange(0, 225, 45))
        ax.tick_params(axis='both', which='major', labelsize=24)
        ax.set_xlim([-180, 180])
        ax.set_ylim([0, 180])

        # Adding top conformation 
        ax.text(-20.0, 180.0+3.0, '\u2074C\u2081', fontsize=30, fontweight='bold', color='black')
        ax.text(-20.0, -20.0, '\u00b9C\u2084', fontsize=30, fontweight='bold', color='black')
       
    
        # Add grid lines
        #ax.xaxis.grid(True)
        #ax.yaxis.grid(True)
        ax.grid(which='major', linewidth=0.8, color='black', linestyle='--', alpha = 0.8)
    
        # Add titles
        if ax == axn.flat[0]:
            ax.set_title('Canonical RUs, vacuum\n', fontweight="bold", fontsize=40)
        elif ax == axn.flat[1]:
            ax.set_title('Non-canonical RUs, vacuum\n', fontweight="bold", fontsize=40)
        elif ax == axn.flat[2]:
            ax.set_title('All RUs, vacuum\n', fontweight="bold", fontsize=40)
        elif ax == axn.flat[3]:
            ax.set_title('\nCanonical RUs, water\n', fontweight="bold", fontsize=40)
        elif ax == axn.flat[4]:
            ax.set_title('\nNon-canonical RUs, water\n', fontweight="bold", fontsize=40)
        elif ax == axn.flat[5]:
            ax.set_title('\nAll RUs, water\n', fontweight="bold", fontsize=40)
        
        # Set legend
        if ax == axn.flat[5]:
            grey_square = mlines.Line2D([], [], color='grey', marker='s', linestyle='None',
                          markersize=16, label='α')
            green_triangle = mlines.Line2D([], [], color='green', marker='^', linestyle='None',
                             markersize=16, label='β')
            ax.legend(handles=[grey_square, green_triangle], fontsize=40, bbox_to_anchor=(1.2, -0.2), loc="lower center", ncol=2)
                
    plt.tight_layout(h_pad=1)
    plt.savefig('Scatter_plot_'+str(TC)+'.pdf')
    
    # Generate freq histogram for the puckering amplitude
    # Vacuum
    a = db_vacuum_canonical_alpha['Q(Total_amplitude)_Cremer&Pople'].to_numpy()
    b = db_vacuum_canonical_beta['Q(Total_amplitude)_Cremer&Pople'].to_numpy()
    Q_db_vacuum_canonical = pd.DataFrame({'Q_alpha': a,
                                          'Q_beta': b})

    a = db_vacuum_non_canonical_alpha['Q(Total_amplitude)_Cremer&Pople'].to_numpy()
    b = db_vacuum_non_canonical_beta['Q(Total_amplitude)_Cremer&Pople'].to_numpy()
    Q_db_vacuum_non_canonical = pd.DataFrame({'Q_alpha': a,
                                              'Q_beta': b})

    Q_db_vacuum_all = pd.DataFrame({'Q_alpha': db_vacuum_alpha['Q(Total_amplitude)_Cremer&Pople'].to_numpy(),
                                    'Q_beta': db_vacuum_beta['Q(Total_amplitude)_Cremer&Pople'].to_numpy()})

    # Water
    a = db_water_canonical_alpha['Q(Total_amplitude)_Cremer&Pople'].to_numpy()
    b = db_water_canonical_beta['Q(Total_amplitude)_Cremer&Pople'].to_numpy()
    Q_db_water_canonical = pd.DataFrame({'Q_alpha': a,
                                         'Q_beta': b})

    a = db_water_non_canonical_alpha['Q(Total_amplitude)_Cremer&Pople'].to_numpy()
    b = db_water_non_canonical_beta['Q(Total_amplitude)_Cremer&Pople'].to_numpy()
    Q_db_water_non_canonical = pd.DataFrame({'Q_alpha': a,
                                             'Q_beta': b})

    Q_db_water_all = pd.DataFrame({'Q_alpha': db_water_alpha['Q(Total_amplitude)_Cremer&Pople'].to_numpy(),
                                   'Q_beta': db_water_beta['Q(Total_amplitude)_Cremer&Pople'].to_numpy()})
    
    counts, bins = np.histogram(Q_db_vacuum_canonical['Q_alpha']) 
    
    # Generate density prob histograms for the Q total amplitude
    fig, axn = plt.subplots(2, 3, sharex=False, sharey=False, figsize=(18,10)) 
    for i, (ax, dataframe, env, type_ru) in enumerate([(axn.flat[0], Q_db_vacuum_canonical, "vacuum", "canonical"), 
                                                       (axn.flat[1], Q_db_vacuum_non_canonical, "vacuum", "non_canonical"), 
                                                       (axn.flat[2], Q_db_vacuum_all, "vacuum", "all"),           
                                                       (axn.flat[3], Q_db_water_canonical, "water", "canonical"),
                                                       (axn.flat[4], Q_db_water_non_canonical, "water", "non_canonical"),          
                                                       (axn.flat[5], Q_db_water_all, "water", "all")]):
    
        # Generate the density plots
        counts_alpha, bins_alpha = np.histogram(dataframe['Q_alpha'])
        prob_alpha = (counts_alpha/sum(counts_alpha))*100
        counts_beta, bins_beta = np.histogram(dataframe['Q_beta'])
        prob_beta = (counts_beta/sum(counts_alpha))*100
        bar_width = 0.4
    
        #ax.bar(bins_alpha[1:], prob_alpha, width = bar_width, align="edge", color = "grey", 
          #     edgecolor = 'black')
        #ax.bar(bins_beta[1:], prob_beta, width = bar_width, align="edge", color = "green", 
         #      edgecolor = 'black')
    
        if ax == axn.flat[5]:
            hist = ax.hist([dataframe['Q_alpha'], dataframe['Q_beta']],
                            bins = 'auto', stacked=False, density=False,
                            color = ['grey','g'], edgecolor = 'black', label=['alpha', 'beta'])
        else:
            hist = ax.hist([dataframe['Q_alpha'], dataframe['Q_beta']],
                            bins = 'auto', stacked=False, density=False,
                            color = ['grey','g'], edgecolor = 'black')
    
        #def adjust_y_axis(x, pos):
            #return (x / (len(mydata) * 1.0)*100
    
        # Set ticks sizes
        ax.tick_params(labelsize=18)
    
        ## Set labels titles
        if ax == axn.flat[0]:
            ax.set_ylabel('Count', fontsize = 18, fontweight='bold') 
        elif ax == axn.flat[3]:
            ax.set_xlabel("Puckering amplitude Q(Å)", fontsize = 18, fontweight='bold')
            ax.set_ylabel('Count', fontsize = 18, fontweight='bold') 
        elif ax == axn.flat[4] or ax == axn.flat[5]:  
            ax.set_xlabel("Puckering amplitude Q(Å)", fontsize = 18, fontweight='bold')
            
        # Add labels for bars
        for container in ax.containers:
            ax.bar_label(container, fmt='%.1f', fontsize=12)
        
        # Set title
        # Add the title to rose diagram
        if ax == axn.flat[0]:
            ax.set_title('Canonical, vacuum', fontweight="bold", fontsize=18)
        elif ax == axn.flat[1]:
            ax.set_title('Non_canonical, vacuum', fontweight="bold", fontsize=18)
        elif ax == axn.flat[2]:
            ax.set_title('All RUs, vacuum', fontweight="bold", fontsize=18)
        elif ax == axn.flat[3]:
            ax.set_title('\nCanonical, water', fontweight="bold", fontsize=18)
        elif ax == axn.flat[4]:
            ax.set_title('\nNon canonical, water', fontweight="bold", fontsize=18)
        elif ax == axn.flat[5]:
            ax.set_title('\nAll RUs, water', fontweight="bold", fontsize=18)
            
        # Set legend
        if ax == axn.flat[5]:
            ax.legend(fontsize=18, bbox_to_anchor=(1.65, 0.00), loc="lower right")
    
    plt.tight_layout()
    plt.savefig('Q_freq_hist_'+str(TC)+'.pdf')     
       
END_SCRIPT

# Run script.py with python3
python3 gen_ring_puckering_csv_graph.py

cd ../ # LEAVE CURRENT TC FOLDER
done








