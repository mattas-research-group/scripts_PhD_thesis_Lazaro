#!/bin/bash

####################################################################
### DESCRIPTION:                                                 ###
### This script will examine each scan log file. Count number    ###
### of structures and if  = 7 then it completed otherwise not    ###
### completed. In case of completion obtain minimum and generate ###
### final gaussian opt file. Generate list of angles and         ###
### relative energies in kcal/mol.                               ###
### In case of not completed output if was not archived or did   ###
### not finished normally in different (.txt) files              ###
####################################################################

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a RU_names=("adenine" "guanine" "cytosine" "thymine" "uracil" "TARC_cbond" "TARC_nbond" "BA_cbond" "BA_nbond" "CA" "melamine")
declare -a RU_ac_names=("adenine_ac" "guanine_ac" "cytosine_ac" "thymine_ac" "uracil_ac" "TARC_cbond_ac" "TARC_nbond_ac" "BA_cbond_ac" "BA_nbond_ac" "CA_ac" "melamine_ac")
declare -a IL_names=("phosphate" "arsenate")

# Create folder for gaussian opt + freq of nucleotides for a+b pathway
mkdir $MD_RUN/nucleotides_a+b/final_gaussian_opt



function completion_status {

# Assign variables names to entries of function
local environment=$1
local log_file=$2
local file_name=$3
local type_sugar=$4
local sugar_ring=$5
local RU=$6
local IL=$7

# Delete previous dat files
rm *.dat

# Find all instances of convergence message and print how many are. If they are 7 then continue otherwise generate a new file to rerun the calculation
grep -o 'Stationary point found' $log_file | wc -l > number_conformers_that_converged.txt

# Assign variable to this file
num_conf_converged=`cat number_conformers_that_converged.txt`

if [ "$num_conf_converged" == 7 ]; then

# Print angles values in different file
for angle_value in $(seq 0 60 360); do
echo $angle_value >> list_torsion_values.dat
done

# Now work on obtain the energies for each scan stationary point and their respective geometries and get the minimum
for structure_from_scan in $(seq 1 1 7); do

# Get line number of each instance of -- Stationary point found.
awk '/-- Stationary point found./{c++} c=='$structure_from_scan'{print NR;exit}' $log_file > line_number_conv_message.txt
line_number_with_conv_message=`cat line_number_conv_message.txt`


# Print line number for next convergence message -- Stationary point found.
if [ "$structure_from_scan" -eq 1 ]; then

# Copy from beginning to first convergence message -- Stationary point found.
sed -n '1,'$line_number_with_conv_message'p' $log_file > structure#"$structure_from_scan"_file.txt 

# Copy energy value for this structure from scan 
#Print last instance of SCF Done:     E(RB3LYP)
awk '/SCF Done/ { ln = FNR } END { print ln }' structure#"$structure_from_scan"_file.txt > line_number_with_energy_structure#"$structure_from_scan".txt
energy_line_number=`cat line_number_with_energy_structure#"$structure_from_scan".txt`

# Copy this line in separated file
sed -n ''$energy_line_number','$energy_line_number'p' structure#"$structure_from_scan"_file.txt > line_with_energy_structure#"$structure_from_scan".txt

# Clean up this file to get only the energy
awk '{ print $5 }' line_with_energy_structure#"$structure_from_scan".txt >> energies_from_scan.txt 


elif [ "$structure_from_scan" -gt 1 ]; then


# Obtain previous structure number
previous_structure_from_scan=$(($structure_from_scan-1))

# Get line number of each instance of -- Stationary point found.
awk '/-- Stationary point found./{c++} c=='$previous_structure_from_scan'{print NR;exit}' $log_file > previous_line_number_conv_message.txt
previous_line_number_with_conv_message=`cat previous_line_number_conv_message.txt`
sed -n ''$previous_line_number_with_conv_message','$line_number_with_conv_message'p' $log_file > structure#"$structure_from_scan"_file.txt 

# Copy energy value for this structure from scan 
#Print last instance of SCF Done:     E(RB3LYP)
awk '/SCF Done/ { ln = FNR } END { print ln }' structure#"$structure_from_scan"_file.txt > line_number_with_energy_structure#"$structure_from_scan".txt
energy_line_number=`cat line_number_with_energy_structure#"$structure_from_scan".txt`

# Copy this line in separated file
sed -n ''$energy_line_number','$energy_line_number'p' structure#"$structure_from_scan"_file.txt > line_with_energy_structure#"$structure_from_scan".txt

# Clean up this file to get only the energy
awk '{ print $5 }' line_with_energy_structure#"$structure_from_scan".txt >> energies_from_scan.txt 

fi

# Print structure number in separated file
echo "$structure_from_scan" >> structure_number_from_scan.txt

done

# Convert energies to kcal/mol
awk '{ printf "%.8e\n", $1*627.503 }' energies_from_scan.txt > energies_from_scan_kcal_mol.txt

# Obtain min from this file
sort -k1 -n energies_from_scan_kcal_mol.txt > energies_from_scan_kcal_mol_sorted.txt
awk 'FNR == 1 {print $1}' energies_from_scan_kcal_mol_sorted.txt > min_energy_kcal_mol.txt
min_energy=`cat min_energy_kcal_mol.txt`

# Repeat min energy
for i in `seq $num_conf_converged`; do
echo "$min_energy" >> min_energy_kcal_mol_repeated.txt
done

# Put all energies and min energy in same file
paste energies_from_scan_kcal_mol.txt min_energy_kcal_mol_repeated.txt > to_calc_rel_energies.txt

# Obtain relative energies in kcal/mol
awk '{ printf "%.8e\n", $1=$1-$2 }' to_calc_rel_energies.txt > list_relative_energies_scan_kcal_mol.dat

# Put file with structures number list and energies list in one
paste structure_number_from_scan.txt energies_from_scan.txt > list_energies_structures.txt

# Order this file in decresing order of the energies (column 2)
sort -k2 -n list_energies_structures.txt > list_sorted_energies_structures.txt

# Print first column and first line of this column to obtain scan global minimum structure#
awk 'FNR == 1 {print $1}' list_sorted_energies_structures.txt > minimum_structure_number.txt
minimum_structure_number=`cat minimum_structure_number.txt`

# Get cartesian coordinates of atoms for minimum structure
#Get last instance of standard orientation line
awk '/Standard orientation/ { ln = FNR } END { print ln }' structure#"$minimum_structure_number"_file.txt > start_line_mol_structure_structure#"$minimum_structure_number".txt
start_line=`cat start_line_mol_structure_structure#"$minimum_structure_number".txt` 

# Get last instance of Rotational constants line
awk '/Rotational constants/ { ln = FNR } END { print ln }' structure#"$minimum_structure_number"_file.txt > last_line_mol_structure_structure#"$minimum_structure_number".txt
last_line=`cat last_line_mol_structure_structure#"$minimum_structure_number".txt`

# Copy cartesian coordinates info
sed -n ''$start_line','$last_line'p' structure#"$minimum_structure_number"_file.txt > cartesian_coord_min.txt

# Print only the last three columns of cartesian coordinates
awk '{ print $4, $5, $6 }' cartesian_coord_min.txt > tmp.txt && mv tmp.txt cartesian_coord_min.txt

# Delete last two lines of this file
head -n -2 cartesian_coord_min.txt > tmp.txt && mv tmp.txt cartesian_coord_min.txt
# Delete first 5 lines of this file
tail -n +6 cartesian_coord_min.txt > tmp.txt && mv tmp.txt cartesian_coord_min.txt

# Justify columns to the right
sed 's/$/\tX/g' cartesian_coord_min.txt |column -t |sed -r 's/([-+]?[0-9.]+)( +)/\2\1/g; s/^  //; s/X$//' > tmp.txt && mv tmp.txt cartesian_coord_min.txt

# Work now on gjf file to obtain the atoms symbols
# Delete first 8 lines
tail -n +9 *.gjf > atoms_symbols.txt

# Delete last two lines
head -n -2 atoms_symbols.txt > tmp.txt && mv tmp.txt atoms_symbols.txt 

# Print just first column of this file
awk '{ print $1 }' atoms_symbols.txt > tmp.txt && mv tmp.txt atoms_symbols.txt

# Print columns 4, 5, 6 of this file
awk '{ print  $1, "                 " }' atoms_symbols.txt > tmp.txt && mv tmp.txt atoms_symbols.txt

# Put all together 
paste atoms_symbols.txt cartesian_coord_min.txt > nucleotide.gjf

# Add charge and multiplicity
sed -i '1 i\-1 1' nucleotide.gjf

# Add instructions
#sed -i '1 i\-1 1' nucleoside.gjf
echo " " > title_gjf_file.txt
echo "Opt and freq optimization for PES minimum (min#"$minimum_structure_number") of "$type_sugar"_"$sugar_ring"_"$ru"_"$IL" in "$environment" environment" >> title_gjf_file.txt
echo " " >> title_gjf_file.txt

# Paste DFT instructions
cat $MD_RUN/instructions_DFT/instructions_b3lyp.txt title_gjf_file.txt nucleotide.gjf > tmp.txt && mv tmp.txt "$file_name"_min.gjf
echo " " >> "$file_name"_min.gjf

# Replace name of chk file in final gjf file
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$file_name'_min.chk/g' "$file_name"_min.gjf

# Delete last 2 lines of this file
head -n -2 "$file_name"_min.gjf > tmp.txt && mv tmp.txt "$file_name"_min.gjf

# Substitute keywords for the case of water opt
if [ "$environment" == "water" ]; then
sed -i 's/scf=maxcycle=1600/scf=maxcycle=1600 scrf=(solvent=water)/g' "$file_name"_min.gjf
fi

# Replace instructions for DFT calculation
sed -i 's/opt=modredundant/opt=calcall freq/g' "$file_name"_min.gjf
sed -i 's/6-31g(d,p)/6-311++g(d,p)/g' "$file_name"_min.gjf

# Delete first line from file
sed -i 1,1d "$file_name"_min.gjf

# Replace RAM and PROC info
sed -i 's/%nprocshared=64/%nprocshared=16/g' "$file_name"_min.gjf
sed -i 's/%mem=160GB/%mem=16GB/g' "$file_name"_min.gjf

# Add empty line at the end
echo "" >> "$file_name"_min.gjf

# Copy file job1.sh script to corresponding folder
cp "$file_name"_min.gjf $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$IL/$RU


if [ "$sugar_ring" == "none" ]; then

# Copy file to corresponding folder
cp "$file_name"_min.gjf $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$IL/$RU

# Copy job1 bash to this folder
cp $MD_RUN/job1.sh $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$IL/$RU

# Replace name of file to calculate in job1
sed -i 's/g16 2deoxy_ribofuranose_2endo_alpha_guanine_phosphate_scan2.gjf/g16 '$file_name'_min.gjf/g' $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$IL/$RU/job1.sh

else

# Copy file to corresponding folder
cp "$file_name"_min.gjf $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$sugar_ring/$IL/$RU

# Copy job1 bash to this folder
cp $MD_RUN/job1.sh $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$sugar_ring/$IL/$RU

# Replace name of file to calculate in job1
sed -i 's/g16 2deoxy_ribofuranose_2endo_alpha_guanine_phosphate_scan2.gjf/g16 '$file_name'_min.gjf/g' $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$sugar_ring/$IL/$RU/job1.sh

fi


rm *.txt






else # Case in which the scan did not finished normally







# Print angles values in different file
for angle_value in $(seq 0 60 360); do
echo $angle_value >> list_torsion_values.dat
done

# Now work on obtain the energies for each scan stationary point and their respective geometries and get the minimum
for structure_from_scan in $(seq 1 1 $num_conf_converged); do

# Get line number of each instance of -- Stationary point found.
awk '/-- Stationary point found./{c++} c=='$structure_from_scan'{print NR;exit}' $log_file > line_number_conv_message.txt
line_number_with_conv_message=`cat line_number_conv_message.txt`


# Print line number for next convergence message -- Stationary point found.
if [ "$structure_from_scan" -eq 1 ]; then

# Copy from beginning to first convergence message -- Stationary point found.
sed -n '1,'$line_number_with_conv_message'p' $log_file > structure#"$structure_from_scan"_file.txt 

# Copy energy value for this structure from scan 
#Print last instance of SCF Done:     E(RB3LYP)
awk '/SCF Done/ { ln = FNR } END { print ln }' structure#"$structure_from_scan"_file.txt > line_number_with_energy_structure#"$structure_from_scan".txt
energy_line_number=`cat line_number_with_energy_structure#"$structure_from_scan".txt`

# Copy this line in separated file
sed -n ''$energy_line_number','$energy_line_number'p' structure#"$structure_from_scan"_file.txt > line_with_energy_structure#"$structure_from_scan".txt

# Clean up this file to get only the energy
awk '{ print $5 }' line_with_energy_structure#"$structure_from_scan".txt >> energies_from_scan.txt 


elif [ "$structure_from_scan" -gt 1 ]; then


# Obtain previous structure number
previous_structure_from_scan=$(($structure_from_scan-1))

# Get line number of each instance of -- Stationary point found.
awk '/-- Stationary point found./{c++} c=='$previous_structure_from_scan'{print NR;exit}' $log_file > previous_line_number_conv_message.txt
previous_line_number_with_conv_message=`cat previous_line_number_conv_message.txt`
sed -n ''$previous_line_number_with_conv_message','$line_number_with_conv_message'p' $log_file > structure#"$structure_from_scan"_file.txt 

# Copy energy value for this structure from scan 
#Print last instance of SCF Done:     E(RB3LYP)
awk '/SCF Done/ { ln = FNR } END { print ln }' structure#"$structure_from_scan"_file.txt > line_number_with_energy_structure#"$structure_from_scan".txt
energy_line_number=`cat line_number_with_energy_structure#"$structure_from_scan".txt`

# Copy this line in separated file
sed -n ''$energy_line_number','$energy_line_number'p' structure#"$structure_from_scan"_file.txt > line_with_energy_structure#"$structure_from_scan".txt

# Clean up this file to get only the energy
awk '{ print $5 }' line_with_energy_structure#"$structure_from_scan".txt >> energies_from_scan.txt 

fi

# Print structure number in separated file
echo "$structure_from_scan" >> structure_number_from_scan.txt

done

# Convert energies to kcal/mol
awk '{ printf "%.8e\n", $1*627.503 }' energies_from_scan.txt > energies_from_scan_kcal_mol.txt

# Obtain min from this file
sort -k1 -n energies_from_scan_kcal_mol.txt > energies_from_scan_kcal_mol_sorted.txt
awk 'FNR == 1 {print $1}' energies_from_scan_kcal_mol_sorted.txt > min_energy_kcal_mol.txt
min_energy=`cat min_energy_kcal_mol.txt`

# Repeat min energy
for i in `seq $num_conf_converged`; do
echo "$min_energy" >> min_energy_kcal_mol_repeated.txt
done

# Put all energies and min energy in same file
paste energies_from_scan_kcal_mol.txt min_energy_kcal_mol_repeated.txt > to_calc_rel_energies.txt

# Obtain relative energies in kcal/mol
awk '{ printf "%.8e\n", $1=$1-$2 }' to_calc_rel_energies.txt > list_relative_energies_scan_kcal_mol.dat

# Put file with structures number list and energies list in one
paste structure_number_from_scan.txt energies_from_scan.txt > list_energies_structures.txt

# Order this file in decresing order of the energies (column 2)
sort -k2 -n list_energies_structures.txt > list_sorted_energies_structures.txt

# Print first column and first line of this column to obtain scan global minimum structure#
awk 'FNR == 1 {print $1}' list_sorted_energies_structures.txt > minimum_structure_number.txt
minimum_structure_number=`cat minimum_structure_number.txt`

# Get cartesian coordinates of atoms for minimum structure
#Get last instance of standard orientation line
awk '/Standard orientation/ { ln = FNR } END { print ln }' structure#"$minimum_structure_number"_file.txt > start_line_mol_structure_structure#"$minimum_structure_number".txt
start_line=`cat start_line_mol_structure_structure#"$minimum_structure_number".txt` 

# Get last instance of Rotational constants line
awk '/Rotational constants/ { ln = FNR } END { print ln }' structure#"$minimum_structure_number"_file.txt > last_line_mol_structure_structure#"$minimum_structure_number".txt
last_line=`cat last_line_mol_structure_structure#"$minimum_structure_number".txt`

# Copy cartesian coordinates info
sed -n ''$start_line','$last_line'p' structure#"$minimum_structure_number"_file.txt > cartesian_coord_min.txt

# Print only the last three columns of cartesian coordinates
awk '{ print $4, $5, $6 }' cartesian_coord_min.txt > tmp.txt && mv tmp.txt cartesian_coord_min.txt

# Delete last two lines of this file
head -n -2 cartesian_coord_min.txt > tmp.txt && mv tmp.txt cartesian_coord_min.txt
# Delete first 5 lines of this file
tail -n +6 cartesian_coord_min.txt > tmp.txt && mv tmp.txt cartesian_coord_min.txt

# Justify columns to the right
sed 's/$/\tX/g' cartesian_coord_min.txt |column -t |sed -r 's/([-+]?[0-9.]+)( +)/\2\1/g; s/^  //; s/X$//' > tmp.txt && mv tmp.txt cartesian_coord_min.txt

# Work now on gjf file to obtain the atoms symbols
# Delete first 8 lines
tail -n +9 *.gjf > atoms_symbols.txt

# Delete last two lines
head -n -2 atoms_symbols.txt > tmp.txt && mv tmp.txt atoms_symbols.txt 

# Print just first column of this file
awk '{ print $1 }' atoms_symbols.txt > tmp.txt && mv tmp.txt atoms_symbols.txt

# Print columns 4, 5, 6 of this file
awk '{ print  $1, "                 " }' atoms_symbols.txt > tmp.txt && mv tmp.txt atoms_symbols.txt

# Put all together 
paste atoms_symbols.txt cartesian_coord_min.txt > nucleotide.gjf

# Add charge and multiplicity
sed -i '1 i\-1 1' nucleotide.gjf

# Add instructions
#sed -i '1 i\-1 1' nucleoside.gjf
echo " " > title_gjf_file.txt
echo "Opt and freq optimization for PES minimum (min#"$minimum_structure_number") of "$type_sugar"_"$sugar_ring"_"$ru"_"$IL" in "$environment" environment" >> title_gjf_file.txt
echo " " >> title_gjf_file.txt

# Paste DFT instructions
cat $MD_RUN/instructions_DFT/instructions_b3lyp.txt title_gjf_file.txt nucleotide.gjf > tmp.txt && mv tmp.txt "$file_name"_min.gjf
echo " " >> "$file_name"_min.gjf

# Replace name of chk file in final gjf file
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$file_name'_min.chk/g' "$file_name"_min.gjf

# Delete last 2 lines of this file
head -n -2 "$file_name"_min.gjf > tmp.txt && mv tmp.txt "$file_name"_min.gjf

# Substitute keywords for the case of water opt
if [ "$environment" == "water" ]; then
sed -i 's/scf=maxcycle=1600/scf=maxcycle=1600 scrf=(solvent=water)/g' "$file_name"_min.gjf
fi

# Replace instructions for DFT calculation
sed -i 's/opt=modredundant/opt=calcall freq/g' "$file_name"_min.gjf
sed -i 's/6-31g(d,p)/6-311++g(d,p)/g' "$file_name"_min.gjf

# Delete first line from file
sed -i 1,1d "$file_name"_min.gjf

# Replace RAM and PROC info
sed -i 's/%nprocshared=64/%nprocshared=16/g' "$file_name"_min.gjf
sed -i 's/%mem=160GB/%mem=16GB/g' "$file_name"_min.gjf

# Add empty line at the end
echo "" >> "$file_name"_min.gjf

# Copy file job1.sh script to corresponding folder
cp "$file_name"_min.gjf $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$IL/$RU


if [ "$sugar_ring" == "none" ]; then

# Copy file to corresponding folder
cp "$file_name"_min.gjf $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$IL/$RU

# Copy job1 bash to this folder
cp $MD_RUN/job1.sh $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$IL/$RU

# Replace name of file to calculate in job1
sed -i 's/g16 2deoxy_ribofuranose_2endo_alpha_guanine_phosphate_scan2.gjf/g16 '$file_name'_min.gjf/g' $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$IL/$RU/job1.sh

else

# Copy file to corresponding folder
cp "$file_name"_min.gjf $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$sugar_ring/$IL/$RU

# Copy job1 bash to this folder
cp $MD_RUN/job1.sh $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$sugar_ring/$IL/$RU

# Replace name of file to calculate in job1
sed -i 's/g16 2deoxy_ribofuranose_2endo_alpha_guanine_phosphate_scan2.gjf/g16 '$file_name'_min.gjf/g' $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$type_sugar/$sugar_ring/$IL/$RU/job1.sh

fi


rm *.txt







fi

}


##########################################################################################################################
############################################## WORK ON EACH FOLDER #######################################################
##########################################################################################################################

# Open corresponding folder
cd nucleotides_a+b/gaussian_scan

for environment in "${environment[@]}"; do

# Create corresponding folder in final_gaussian_opt folder
mkdir $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment

cd $environment

for TC in "${TC_names[@]}"; do

# Create corresponding folder in final_gaussian_opt folder
mkdir $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$TC

cd $TC

if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

for ring_conf_fura in "${Ring_Conf_furanose[@]}"; do

# Create corresponding folder in final_gaussian_opt folder
mkdir $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$TC/$ring_conf_fura

cd $ring_conf_fura

for IL in "${IL_names[@]}"; do

# Create corresponding folder in final_gaussian_opt folder
mkdir $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$TC/$ring_conf_fura/$IL

cd $IL

for RU in "${RU_names[@]}"; do

# Create corresponding folder in final_gaussian_opt folder
mkdir $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$TC/$ring_conf_fura/$IL/$RU

cd $RU

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, ring conf $ring_conf_fura, IL $IL and RU $RU"

# Delete if there is core file large file for not conclusion of calculation
find . -type f -name core\* -exec rm {} \;

# Define input file for function that calculates the thermo equilibrium
log_file_prefix="$TC"_"$ring_conf_fura"_"$RU"_"$IL"

# Define input file for function that calculates the thermo equilibrium
unique_opt_file="$log_file_prefix"_scan.log
second_opt_file="$log_file_prefix"_scan2.log
third_opt_file="$log_file_prefix"_scan3.log

if [ -f "$third_opt_file" ]; then

log_file=$third_opt_file

echo "$environment/$TC/$ring_conf_fura/$IL/$RU has a 3rd scan file" >> $MD_RUN/nucleosides/final_gaussian_opt/$environment/number_of_scans_"$environment".txt

else

if [ -f "$second_opt_file" ]; then

log_file=$second_opt_file

echo "$environment/$TC/$ring_conf_fura/$IL/$RU has a 2nd scan file" >> $MD_RUN/nucleosides/final_gaussian_opt/$environment/number_of_scans_"$environment".txt 

else

log_file=$unique_opt_file

echo "$environment/$TC/$ring_conf_fura/$IL/$RU has 1 scan file" >> $MD_RUN/nucleosides/final_gaussian_opt/$environment/number_of_scans_"$environment".txt

fi

fi

# Run function that analyses the scan file and generates info and min for final gaussian opt
completion_status $environment $log_file $log_file_prefix $TC $ring_conf_fura $RU $IL


cd ../
done

cd ../
done

cd ../
done

elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then

for ring_conf_pyra in "${Ring_Conf_pyranose[@]}"; do

# Create corresponding folder in final_gaussian_opt folder
mkdir $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$TC/$ring_conf_pyra

cd $ring_conf_pyra

for IL in "${IL_names[@]}"; do

# Create corresponding folder in final_gaussian_opt folder
mkdir $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$TC/$ring_conf_pyra/$IL

cd $IL

for RU in "${RU_names[@]}"; do

# Create corresponding folder in final_gaussian_opt folder
mkdir $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$TC/$ring_conf_pyra/$IL/$RU

cd $RU

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, ring conf $ring_conf_pyra, IL $IL and RU $RU"

# Delete if there is core file large file for not conclusion of calculation
find . -type f -name core\* -exec rm {} \;

# Define input file for function that calculates the thermo equilibrium
log_file_prefix="$TC"_"$ring_conf_pyra"_"$RU"_"$IL"

# Define input file for function that calculates the thermo equilibrium
unique_opt_file="$log_file_prefix"_scan.log
second_opt_file="$log_file_prefix"_scan2.log
third_opt_file="$log_file_prefix"_scan3.log

if [ -f "$third_opt_file" ]; then

log_file=$third_opt_file

echo "$environment/$TC/$ring_conf_pyra/$IL/$RU has a 3rd scan file" >> $MD_RUN/nucleosides/final_gaussian_opt/$environment/number_of_scans_"$environment".txt

else

if [ -f "$second_opt_file" ]; then

log_file=$second_opt_file

echo "$environment/$TC/$ring_conf_pyra/$IL/$RU has a 2nd scan file" >> $MD_RUN/nucleosides/final_gaussian_opt/$environment/number_of_scans_"$environment".txt 

else

log_file=$unique_opt_file

echo "$environment/$TC/$ring_conf_pyra/$IL/$RU has 1 scan file" >> $MD_RUN/nucleosides/final_gaussian_opt/$environment/number_of_scans_"$environment".txt

fi

fi

# Run function that analyses the scan file and generates info and min for final gaussian opt
completion_status $environment $log_file $log_file_prefix $TC $ring_conf_pyra $RU $IL




cd ../
done

cd ../
done

cd ../
done

elif [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then

for IL in "${IL_names[@]}"; do

# Create corresponding folder in final_gaussian_opt folder
mkdir $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$TC/$IL

cd $IL

for RU in "${RU_names[@]}"; do

# Create corresponding folder in final_gaussian_opt folder
mkdir $MD_RUN/nucleotides_a+b/final_gaussian_opt/$environment/$TC/$IL/$RU

cd $RU

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, none ring conf, IL $IL and RU $RU"

ring_conf="none"


# Define input file for function that calculates the thermo equilibrium
log_file_prefix="$TC"_"$RU"_"$IL"

# Define input file for function that calculates the thermo equilibrium
unique_opt_file="$log_file_prefix"_scan.log
second_opt_file="$log_file_prefix"_scan2.log
third_opt_file="$log_file_prefix"_scan3.log

if [ -f "$third_opt_file" ]; then

log_file=$third_opt_file

echo "$environment/$TC/$IL/$RU has a 3rd scan file" >> $MD_RUN/nucleosides/final_gaussian_opt/$environment/number_of_scans_"$environment".txt

else

if [ -f "$second_opt_file" ]; then

log_file=$second_opt_file

echo "$environment/$TC/$IL/$RU has a 2nd scan file" >> $MD_RUN/nucleosides/final_gaussian_opt/$environment/number_of_scans_"$environment".txt 

else

log_file=$unique_opt_file

echo "$environment/$TC/$IL/$RU has 1 scan file" >> $MD_RUN/nucleosides/final_gaussian_opt/$environment/number_of_scans_"$environment".txt

fi

fi

# Run function that analyses the scan file and generates info and min for final gaussian opt
completion_status $environment $log_file $log_file_prefix $TC $ring_conf $RU $IL


# Delete if there is core file large file for not conclusion of calculation
find . -type f -name core\* -exec rm {} \;



cd ../
done

cd ../
done

fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT
























