#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a RU_names=("adenine" "guanine" "cytosine" "thymine" "uracil" "TARC_cbond" "TARC_nbond" "BA_cbond" "BA_nbond" "CA" "melamine")
declare -a IL_names=("phosphate" "arsenate")


#################################################################
##### CREATE FUNCTION THAT ASSEMBLE THE TC AND RU TO CREATE #####
##### THE CORRESPONDING NUCLEOSIDE IN EACH ENVIRONMENT      #####
#################################################################

function create_nucleotide_cd {

# Note1: no_int_coord_base.gzmat and no_int_coord_sugar.gzmat are the files with only the z matrix

# Declare corresponding variables
local environment=$1
local type_sugar=$2
local sugar_ring=$3
local ru=$4
local IL=$5

# Give name to trifunctional connector file
if [ "$type_sugar" == "glycerol" ] || [ "$type_sugar" == "glyceric_acid" ]; then
tc="$type_sugar"_"$IL"
else
tc="$type_sugar"_"$sugar_ring"_"$IL"
fi

# Obtain type of anomer
if [ "$type_sugar" != "glycerol" ] && [ "$type_sugar" != "glyceric_acid" ]; then
cut -d "_" -f2 <<< "$sugar_ring" > type_anomer.txt
type_anomer=`cat type_anomer.txt`
fi

# Convert RU log/mol file to gzmat
if [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ] || [ "$ru" == "cytosine" ] || [ "$ru" == "thymine" ] || [ "$ru" == "uracil" ] || [ "$ru" == "CA" ]; then
obabel $ru.log -O $ru.gzmat
else
obabel $ru.mol -O $ru.gzmat
fi

# Convert TC-IL log file to gzmat
unique_opt_file="$tc"_min.log
second_opt_file="$tc"_min2.log
third_opt_file="$tc"_min3.log

if [ -f "$third_opt_file" ]; then

log_file=$third_opt_file

obabel $log_file.log -O $tc.gzmat

else

if [ -f "$second_opt_file" ]; then

log_file=$second_opt_file

obabel $log_file.log -O $tc.gzmat


else

log_file=$unique_opt_file

obabel $log_file.log -O $tc.gzmat

fi

fi

################ Insert options to authomatically change labels of atoms to delete and hemiacetal O and anomeric carbon for sugar ########################

if [ "$type_sugar" == "2deoxy_ribofuranose" ]; then

anomer=8
hemiacetal_O=4
oxygen=16
hydrogen=17

elif [ "$type_sugar" == "ribofuranose" ]; then

anomer=8
hemiacetal_O=4
oxygen=18
hydrogen=19

elif [ "$type_sugar" == "threose" ]; then

anomer=5
hemiacetal_O=2
oxygen=14
hydrogen=15

elif [ "$type_sugar" == "2deoxy_ribopyranose" ]; then

anomer=1
hemiacetal_O=2
oxygen=16
hydrogen=17

elif [ "$type_sugar" == "ribopyranose" ]; then

anomer=1
hemiacetal_O=2
oxygen=18
hydrogen=19

elif [ "$type_sugar" == "glycerol" ]; then

anomer=1
hemiacetal_O=2
oxygen=12
hydrogen=13

elif [ "$type_sugar" == "glyceric_acid" ]; then
anomer=9
hemiacetal_O=1
oxygen=11
hydrogen=12

fi


#################################################################
####################### WORK ON TC-IL FILE ######################
#################################################################

# Declare labels of anomeric carbon
#read -p "insert label of anomeric carbon as appears in gaussview:" anomer ####### ************************(THIS LINE CAN CHANGE DEPENDING ON THE TC)***********************

# Declare label of hemiacetalic oxygen for beta or 
#read -p "insert label of hemiacetalic oxygen as appears in gaussview:" hemiacetal_O ####### ************************(THIS LINE CAN CHANGE DEPENDING ON THE TC)***********************

if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ]; then
atom_before_hemiacetal=$(($hemiacetal_O-1))
atom_before_before_hemiacetal=$(($hemiacetal_O-2))
elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] || [ "$type_sugar" == "glycerol" ] || [ "$type_sugar" == "glyceric_acid" ]; then
atom_before_hemiacetal=$(($hemiacetal_O+1))
atom_before_before_hemiacetal=$(($hemiacetal_O+2))
fi

# First declare the labels of the atoms to be deleted from the TC (in sugar case should be an oxygen and hydrogen)
#read -p "insert labels of oxygen and hydrogen (OH) to be deleted as appears in gaussview(e.g. 17 18):" oxygen hydrogen ########## ************************(THIS LINE CAN CHANGE DEPENDING ON THE TC)***********************

# Obtain the lines number for these atoms
line_oxygen=$(($oxygen+6))
line_hydrogen=$(($hydrogen+6))

# Delete lines in matrix that correspond to the O and H in anomeric position
sed ''$line_oxygen','$line_hydrogen'd' $tc.gzmat > sugar_sin_OH.gzmat



########### Calculate total number of atoms in sugar ############

# Delete first headings of this file
sed '1,6d' $tc.gzmat > sinheadings_sugar.gzmat

# Delete values the internal coordinates
sed -i '1,/Variables:/!d' sinheadings_sugar.gzmat
sed -i '$d' sinheadings_sugar.gzmat

# Count number of lines in this file
< sinheadings_sugar.gzmat wc -l > number_atoms_sugar.txt

atoms_sugar=`cat number_atoms_sugar.txt`
atoms_sugar_ion=$(($atoms_sugar-2))
lines_plus_atoms_in_sugar_ion=$(($atoms_sugar_ion+6))
starting_line_int_coord=$(($atoms_sugar_ion+7))

################################################################

# Copy lines with values of internal coordinates to another file
sed -n ''$starting_line_int_coord',$p' sugar_sin_OH.gzmat > internal_coordinates_values_sugar_sinOH.txt 
# Eliminate last empty line
sed -i '$ d' internal_coordinates_values_sugar_sinOH.txt    

#### Delete lines with values of the distances, angles and dihedrals for the atoms deleted in previous step ###
sed -i '/r'$oxygen'/d' internal_coordinates_values_sugar_sinOH.txt
sed -i '/a'$oxygen'/d' internal_coordinates_values_sugar_sinOH.txt
sed -i '/d'$oxygen'/d' internal_coordinates_values_sugar_sinOH.txt
sed -i '/r'$hydrogen'/d' internal_coordinates_values_sugar_sinOH.txt
sed -i '/a'$hydrogen'/d' internal_coordinates_values_sugar_sinOH.txt
sed -i '/d'$hydrogen'/d' internal_coordinates_values_sugar_sinOH.txt
#grep -i "a'$hydrogen'" internal_coordinates_values_sugar_sinOH.txt  


############# Work on correcting the rlabel, alabel and dlabel in the file without internal coordinates values ################
# Delete values the internal coordinates
sed -i '1,/Variables:/!d' sugar_sin_OH.gzmat

# Delete last line of this file 
sed -i '$d' sugar_sin_OH.gzmat

for lines in $(seq 7 1 $lines_plus_atoms_in_sugar_ion)

do

if (("$lines" >= "$line_oxygen")); then

# Obtain old labels
old_label=$(($lines-4))
new_label=$(($lines-6))
# Substitute old lines by new ones
sed -i ''$lines's/r'$old_label'/r'$new_label'/' sugar_sin_OH.gzmat
sed -i ''$lines's/a'$old_label'/a'$new_label'/' sugar_sin_OH.gzmat
sed -i ''$lines's/d'$old_label'/d'$new_label'/' sugar_sin_OH.gzmat

fi
done

############# Work on correcting the rlabel, alabel and dlabel in the file with internal coordinates values ################
for labels in $(seq 2 1 $atoms_sugar_ion)

do

echo "r$labels=" >> labels_sugar_fix.txt
echo "a$labels=" >> labels_sugar_fix.txt
echo "d$labels=" >> labels_sugar_fix.txt

done

# Eliminate innecesary lines from the previous file
sed -i '/a2/d' labels_sugar_fix.txt
sed -i '/d2/d' labels_sugar_fix.txt
sed -i '/d3/d' labels_sugar_fix.txt

# Remove heading
sed -i '1,1d' internal_coordinates_values_sugar_sinOH.txt

# Eliminate everything before the = sign in old int coord file from the base and eliminate first line
cut -f2 -d"=" internal_coordinates_values_sugar_sinOH.txt > internal_coordinates_values_sugar_bef_nucleo.txt

# Paste new labels and values for internal coordinates
paste labels_sugar_fix.txt internal_coordinates_values_sugar_bef_nucleo.txt > new_internal_coordinates_sugar_bef_nucleo.txt

# Reduce spaces between two columns
sed -i $'s/\t//g' new_internal_coordinates_sugar_bef_nucleo.txt

tr -s " " < new_internal_coordinates_sugar_bef_nucleo.txt > internal_coordinates_sugar_sinOH_before.gzmat
# Delete previous file
rm internal_coordinates_values_sugar_sinOH.txt

# Add heading to internal coordinates values
sed -i '1 i Variables:' internal_coordinates_sugar_sinOH_before.gzmat

cp internal_coordinates_sugar_sinOH_before.gzmat internal_coordinates_values_sugar_sinOH.txt ### THIS IS THE FINAL FILE FOR INTERNAL COORDINATES VALUES OF THE SUGAR

# Change the names
mv sugar_sin_OH.gzmat no_int_coord_sugar_sinOH.gzmat ### THIS IS THE FINAL FILE FOR THE Z MATRIX OF SUGAR WITHOUT THE INT COORD VALUES

rm sinheadings_sugar.gzmat number_atoms_sugar.txt 


#################################################################
####################### WORK ON RUs FILE ########################
#################################################################





}


############################################## WORK ON EACH FOLDER #######################################################

# Open corresponding folder
cd nucleotides_c+d/gaussian_scan

for environment in "${environment[@]}"; do

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
echo "working on $environment, TC $TC, ring conf $ring_conf_fura, IL $IL and RU $RU"

# Run function to create the nucleotides for c+d path
create_nucleotide_cd $environment $TC $ring_conf_fura $IL $RU 






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
echo "working on $environment, TC $TC, ring conf $ring_conf_pyra, IL $IL and RU $RU"

# Time to copy the corresponding files from the nucleosides folders and the IL folders
create_nucleotide_cd $environment $TC $ring_conf_pyra $IL $RU 







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
echo "working on $environment, TC $TC, none ring conf, IL $IL and RU $RU"

ring_conf="none"

# Time to copy the corresponding files from the nucleosides folders and the IL folders
create_nucleotide_cd $environment $TC $ring_conf $IL $RU





cd ../
done

cd ../
done








fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT
