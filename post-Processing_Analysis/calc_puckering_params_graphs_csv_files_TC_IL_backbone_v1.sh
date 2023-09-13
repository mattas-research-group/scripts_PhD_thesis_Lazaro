#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a IL_names=("phosphate" "arsenate")


###################################################################
##                                                               ##
## Function to obtain the generalized puckering parameters from  ##
## Cremer & Pople  (Cremer&Pople (1975)) for 5- and 6-member     ##
## rings.                                                        ##                                                 
##                                                               ##
###################################################################

function calc_puckering_param {

# Load function input variables
local TC=$1

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
    if theta < 0:
        theta = 360 + theta
        
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
        summary_fura_ring_conf_from_C_P.write("Twist(3T2)")
    elif phase_angle > 0 and phase_angle < 18:
        summary_fura_ring_conf_from_C_P.write("3T2---3E(C3'-endo)(North)")
    elif phase_angle == 18:
        summary_fura_ring_conf_from_C_P.swrite("3E(C3'-endo)(North)")
    elif phase_angle > 18 and phase_angle < 36:
        summary_fura_ring_conf_from_C_P.write("3E---3T4")
    elif phase_angle == 36:
        summary_fura_ring_conf_from_C_P.swrite("3T4")
    elif phase_angle > 36 and phase_angle < 54:
        summary_fura_ring_conf_from_C_P.write("3T4---E4(C4'-exo)")
    elif phase_angle == 54:
        summary_fura_ring_conf_from_C_P.swrite("E4(C4'-exo)")    
    elif phase_angle > 54 and phase_angle < 72:
        summary_fura_ring_conf_from_C_P.write("E4(C4'-exo)---0T4")
    elif phase_angle == 72:
        summary_fura_ring_conf_from_C_P.swrite("0T4")        
    elif phase_angle > 72 and phase_angle < 90:
        summary_fura_ring_conf_from_C_P.write("0T4---0E(O4'-endo)")
    elif phase_angle == 90:
        summary_fura_ring_conf_from_C_P.swrite("0E(O4'-endo)(East)")    
    elif phase_angle > 90 and phase_angle < 108:
        summary_fura_ring_conf_from_C_P.write("0E---0T1")
    elif phase_angle == 108:
        summary_fura_ring_conf_from_C_P.write("0T1")    
    elif phase_angle > 108 and phase_angle < 126:
        summary_fura_ring_conf_from_C_P.write("0T1---E1(C1'-exo)")
    elif phase_angle == 126:
        summary_fura_ring_conf_from_C_P.write("E1(C1'-exo)")       
    elif phase_angle > 126 and phase_angle < 144:
        summary_fura_ring_conf_from_C_P.write("E1---2T1")
    elif phase_angle == 144:
        summary_fura_ring_conf_from_C_P.write("2T1")   
    elif phase_angle > 144 and phase_angle < 162:
        summary_fura_ring_conf_from_C_P.write("2T1---2E(C2'-endo)")
    elif phase_angle == 162:
        summary_fura_ring_conf_from_C_P.write("2E(C2'-endo)(South)")         
    elif phase_angle > 162 and phase_angle < 180:
        summary_fura_ring_conf_from_C_P.write("2E---2T3(South)")
    elif phase_angle == 180:
        summary_fura_ring_conf_from_C_P.write("2T3(South)")
    elif phase_angle > 180 and phase_angle < 198:
        summary_fura_ring_conf_from_C_P.write("2T3---E3(C3'-exo)(South)")
    elif phase_angle == 198:
        summary_fura_ring_conf_from_C_P.write("E3(C3'-exo)(South)")
    elif phase_angle > 198 and phase_angle < 216:
        summary_fura_ring_conf_from_C_P.write("E3---4T3")
    elif phase_angle == 216:
        summary_fura_ring_conf_from_C_P.write("4T3")
    elif phase_angle > 216 and phase_angle < 234:
        summary_fura_ring_conf_from_C_P.write("4T3---4E(C4'-endo)")
    elif phase_angle == 234:
        summary_fura_ring_conf_from_C_P.write("4E(C4'-endo)")
    elif phase_angle > 234 and phase_angle < 252:
        summary_fura_ring_conf_from_C_P.write("4E(C4'-endo)---4T0")
    elif phase_angle == 252:
        summary_fura_ring_conf_from_C_P.write("4T0")
    elif phase_angle > 252 and phase_angle < 270:
        summary_fura_ring_conf_from_C_P.write("4T0---E0(C4'-exo)")
    elif phase_angle == 270:
        summary_fura_ring_conf_from_C_P.write("E0(C4'-exo)(West)")  
    elif phase_angle > 270 and phase_angle < 288:
        summary_fura_ring_conf_from_C_P.write("E0---1T0")
    elif phase_angle == 288:
        summary_fura_ring_conf_from_C_P.write("1T0")  
    elif phase_angle > 288 and phase_angle < 306:
        summary_fura_ring_conf_from_C_P.write("1T0---1E(C1'-endo)")
    elif phase_angle == 306:    
        summary_fura_ring_conf_from_C_P.write("1E(C1'-endo)")
    elif phase_angle > 306 and phase_angle < 324: 
       summary_fura_ring_conf_from_C_P.write("1E(C1'-endo)---1T2")
    elif phase_angle == 324:    
       summary_fura_ring_conf_from_C_P.write("1T2")
    elif phase_angle > 324 and phase_angle < 342:
       summary_fura_ring_conf_from_C_P.write("1T2--E2(C2'-exo)")
    elif phase_angle == 342:    
       summary_fura_ring_conf_from_C_P.write("E2(C2'-exo)(North)")
    elif phase_angle > 342 and phase_angle < 360:
       summary_fura_ring_conf_from_C_P.write("E2(C2'-exo)---3T2")


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
local IL=$4

# Obtain cartesian coordinates file
log_file=`cat final_list_gaussian.txt`

# Convert log file to gjf to obtain XYZ coord
obabel $log_file -O cartesian_coordinates.gjf

# Delete initial lines, charge and multiplicity to leave only the cartesian coordinates and atoms
sed -i '1,6d' cartesian_coordinates.gjf

# Obtain only the x, y, z coordinates columns
awk '{print $2, $3, $4}' cartesian_coordinates.gjf | column -t > cartesian_coordinates.txt

rm cartesian_coordinates.gjf


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

calc_puckering_param $TC

fi


}








function w_results_to_analysis_folder {

# Load function input variables
local environment=$1
local TC=$2
local ring_conf=$3
local IL=$4

# Send results to general files outside in nucleosides folder
if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

# Load Cremer & Pople parameters
ring_conformation_C_P=`cat 5-member_ring_conformation_C_P.txt`
q2_amplitude_C_P=`cat amplitude_q2_from_C_P.txt`
total_amplitude_Q_C_P=`cat total_amplitude_Q_from_C_P.txt`
phase_angle_phi2_C_P=`cat phase_angle_Phi2.txt`

# Output variables in general files inside the environment folder depending on the ring type
echo "$ring_conformation_C_P" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/final_ring_conformation_C_P_furanoses.txt
echo "$q2_amplitude_C_P" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/q2_C_P_furanoses.txt
echo "$total_amplitude_Q_C_P" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/total_amplitude_Q_C_P_furanoses.txt
echo "$phase_angle_phi2_C_P" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/phase_angle_phi2_Q_C_P_furanoses.txt
echo "$TC" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/TC_names_furanoses.txt
echo "$ring_conf" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/initial_ring_conf_furanoses.txt
echo "$environment" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/environment_list_furanoses.txt
echo "$IL" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/IL_names_furanoses.txt



elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then

# Load C & P parameters
ring_conformation_pyra_C_P=`cat 6-member_ring_conformation_C_P_from_theta.txt`
q2_amplitude_pyra_C_P=`cat amplitude_q2_from_C_P.txt`
q3_amplitude_pyra_C_P=`cat amplitude_q3_from_C_P.txt`
total_amplitude_Q_pyra_C_P=`cat total_amplitude_Q_from_C_P.txt`
phase_angle_phi2_pyra_C_P=`cat phase_angle_Phi2.txt`
phase_angle_theta_C_P=`cat phase_angle_theta.txt`

# Output variables in general files inside the environment folder depending on the ring type
echo "$ring_conformation_pyra_C_P" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/final_ring_conformation_C_P_pyranoses.txt
echo "$q2_amplitude_pyra_C_P" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/q2_C_P_pyranoses.txt
echo "$q3_amplitude_pyra_C_P" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/q3_C_P_pyranoses.txt
echo "$total_amplitude_Q_pyra_C_P" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/total_amplitude_Q_C_P_pyranoses.txt
echo "$phase_angle_phi2_pyra_C_P" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/phase_angle_phi2_Q_C_P_pyranoses.txt
echo "$phase_angle_theta_C_P" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/phase_angle_theta_Q_C_P_pyranoses.txt
echo "$TC" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/TC_names_pyranoses.txt
echo "$ring_conf" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/initial_ring_conf_pyranoses.txt
echo "$environment" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/environment_list_pyranoses.txt
echo "$IL" >> $MD_RUN/post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC/IL_names_pyranoses.txt

fi

}













############################################## WORK ON EACH FOLDER #######################################################

# Create folder for analysis and post-processing
mkdir post_proc_summary_results_csv_graphs/sugar_puckering


# Create IL folder inside the post-processing folder
for IL in "${IL_names[@]}"; do
mkdir post_proc_summary_results_csv_graphs/sugar_puckering/$IL

# Create TC folder inside the post-processing folder
for TC in "${TC_names[@]}"; do
mkdir post_proc_summary_results_csv_graphs/sugar_puckering/$IL/$TC
done
done


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

# Run function to obtain the ring atoms coordinates
obtain_ring_atoms_coord $environment $TC $ring_conf_fura $IL 

# Write results to general files in analysis folder
w_results_to_analysis_folder $environment $TC $ring_conf_fura $IL








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

# Run function to obtain the ring atoms coordinates
obtain_ring_atoms_coord $environment $TC $ring_conf_pyra $IL

# Write results to general files in analysis folder
w_results_to_analysis_folder $environment $TC $ring_conf_pyra $IL







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

cd post_proc_summary_results_csv_graphs/sugar_puckering

for IL in "${IL_names[@]}"; do

cd $IL

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
export IL

cat > gen_ring_puckering_csv_graph.py <<'END_SCRIPT'

# Import environmental variables from python libraries
import sys, os
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from numpy import loadtxt
import matplotlib
import math
from scipy.stats import multivariate_normal
#from astropy.stats import circcorrcoef
#from astropy import units as u
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

# Import environmental variables from bash script for angles in degrees
TC = os.environ.get("TC")
IL = os.environ.get("IL")

# Work on each TC

if TC in ["2deoxy_ribofuranose", "ribofuranose", "threose"]:

    ############## Create csv file for furanose ring data #################
    environment = np.loadtxt("environment_list_furanoses.txt", dtype='str')
    TC_fura = np.loadtxt("TC_names_furanoses.txt", dtype='str')
    initial_ring_conf_fura = np.loadtxt("initial_ring_conf_furanoses.txt", dtype='str')
    final_ring_conf_fura = np.loadtxt("final_ring_conformation_C_P_furanoses.txt", dtype='str')
    q2_fura_C_P = np.loadtxt("q2_C_P_furanoses.txt", dtype='float')
    Q_fura_C_P = np.loadtxt("total_amplitude_Q_C_P_furanoses.txt", dtype='float')
    phi2_fura_C_P = np.loadtxt("phase_angle_phi2_Q_C_P_furanoses.txt", dtype='float')
    
    puckering_df = pd.DataFrame({'Environment': environment, 'TC': TC_fura, 
                                 'Initial_ring_conf': initial_ring_conf_fura, 
                                 'q2(amplitude)_Cremer&Pople': q2_fura_C_P, 
                                 'Q(Total_amplitude)_Cremer&Pople': Q_fura_C_P, 
                                 'phi2(Phase_angle)_Cremer&Pople':phi2_fura_C_P, 
                                 'Final_ring_conf_Cremer&Pople': final_ring_conf_fura})

    # Convert pandas dataframe to csv file
    puckering_df.to_csv('db_'+str(TC)+'_'+str(IL)+'_ring_puckering_params.csv', index=False)

else:

    ############## Create csv file for pyranose ring data #################
    environment = np.loadtxt("environment_list_pyranoses.txt", dtype='str')
    TC_pyra = np.loadtxt("TC_names_pyranoses.txt", dtype='str')
    initial_ring_conf_pyra = np.loadtxt("initial_ring_conf_pyranoses.txt", dtype='str')
    final_ring_conf_pyra = np.loadtxt("final_ring_conformation_C_P_pyranoses.txt", dtype='str')    
    q2_pyra_C_P = np.loadtxt("q2_C_P_pyranoses.txt", dtype='float')
    q3_pyra_C_P = np.loadtxt("q3_C_P_pyranoses.txt", dtype='float')
    Q_pyra_C_P = np.loadtxt("total_amplitude_Q_C_P_pyranoses.txt", dtype='float')
    phi2_pyra_C_P = np.loadtxt("phase_angle_phi2_Q_C_P_pyranoses.txt", dtype='str')
    theta_pyra_C_P = np.loadtxt("phase_angle_theta_Q_C_P_pyranoses.txt", dtype='float')

    puckering_df = pd.DataFrame({'Environment':environment, 'TC': TC_pyra, 
                                 'Initial_ring_conf': initial_ring_conf_pyra, 
                                 'q2(amplitude)_Cremer&Pople': q2_pyra_C_P, 
                                 'q3(amplitude)_Cremer&Pople': q3_pyra_C_P, 
                                 'Q(Total_amplitude)_Cremer&Pople': Q_pyra_C_P, 
                                 'phi2(phase_angle)_Cremer&Pople':phi2_pyra_C_P, 
                                 'theta(phase_angle)_Cremer&Pople':theta_pyra_C_P, 
                                 'Final_ring_conf_Cremer&Pople': final_ring_conf_pyra})
    # Convert pandas dataframe to csv file
    puckering_df.to_csv('db_'+str(TC)+'_'+str(IL)+'_ring_puckering_params.csv', index=False)




END_SCRIPT

# Run script with python3
python3 gen_ring_puckering_csv_graph.py

cd ../ # LEAVE TC FOLDER
done
cd ../ # LEAVE IL FOLDER
done




