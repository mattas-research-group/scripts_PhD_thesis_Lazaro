#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# Insert the names of the files for the TC and the RUs
#read -p "insert name of log file for the Trifunctional Connector (TC) without extension:" tc
read -p "insert name of log file for the Recognition unit (RU) without extension:" ru
#read -p "insert type of anomer in case that it is (e.g. beta, or alpha, or none if there is not anomer)" type_anomer
read -p "insert type of sugar in case that it is (e.g. 2deoxy_ribofuranose, or ribofuranose, or 2deoxy_ribopyranose, or ribopyranose, or threose, or none if the TC is not a sugar)" type_sugar
read -p "insert conformation of the sugar ring (e.g. 2endo_alpha, or 2endo_beta, or D1C4_alpha, or D4C1_beta)" sugar_ring
read -p "input environment (e.g. vacuum or water)" environment

if [ "$sugar_ring" == "none" ]; then

cd $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/$ru

else

cd $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/$sugar_ring/$ru

fi


# Find all instances of convergence message and print how many are. If they are 7 then continue otherwise generate a new file to rerun the calculation
if [ "$sugar_ring" == "none" ]; then

grep -o 'Stationary point found' "$type_sugar"_"$ru"_scan2.log | wc -l > number_conformers_that_converged.txt

else

grep -o 'Stationary point found' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.log | wc -l > number_conformers_that_converged.txt

fi

# Assign variable to this file
num_conf_converged=`cat number_conformers_that_converged.txt`

if [ "$num_conf_converged" == 7 ]; then

#echo "The gaussian scan of the nucleoside of $type_sugar with $sugar_ring sugar ring conformation and RU $ru in $environment environment finished normally" >> $MD_RUN/list_calc_that_finished_normally.log

# Print angles values in different file
#for angle_value in $(seq 0 60 360); do
#echo $angle_value >> list_torsion_values.dat
#done

# Now work on obtain the energies for each scan stationary point and their respective geometries and get the minimum
for structure_from_scan in $(seq 1 1 7); do

# Get line number of each instance of -- Stationary point found.
if [ "$sugar_ring" == "none" ]; then
awk '/-- Stationary point found./{c++} c=='$structure_from_scan'{print NR;exit}' "$type_sugar"_"$ru"_scan2.log > line_number_conv_message.txt
line_number_with_conv_message=`cat line_number_conv_message.txt`
else
awk '/-- Stationary point found./{c++} c=='$structure_from_scan'{print NR;exit}' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.log > line_number_conv_message.txt
line_number_with_conv_message=`cat line_number_conv_message.txt`
fi

# Print line number for next convergence message -- Stationary point found.
if [ "$structure_from_scan" -eq 1 ]; then

# Copy from beginning to first convergence message -- Stationary point found.
if [ "$sugar_ring" == "none" ]; then
sed -n '1,'$line_number_with_conv_message'p' "$type_sugar"_"$ru"_scan2.log > structure#"$structure_from_scan"_file.txt 
else
sed -n '1,'$line_number_with_conv_message'p' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.log > structure#"$structure_from_scan"_file.txt 
fi

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
if [ "$sugar_ring" == "none" ]; then
awk '/-- Stationary point found./{c++} c=='$previous_structure_from_scan'{print NR;exit}' "$type_sugar"_"$ru"_scan2.log > previous_line_number_conv_message.txt
previous_line_number_with_conv_message=`cat previous_line_number_conv_message.txt`
sed -n ''$previous_line_number_with_conv_message','$line_number_with_conv_message'p' "$type_sugar"_"$ru"_scan2.log > structure#"$structure_from_scan"_file.txt 
else
awk '/-- Stationary point found./{c++} c=='$previous_structure_from_scan'{print NR;exit}' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.log > previous_line_number_conv_message.txt
previous_line_number_with_conv_message=`cat previous_line_number_conv_message.txt`
sed -n ''$previous_line_number_with_conv_message','$line_number_with_conv_message'p' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.log > structure#"$structure_from_scan"_file.txt 
fi

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
if [ "$sugar_ring" == "none" ]; then
tail -n +9 "$type_sugar"_"$ru"_scan.gjf > atoms_symbols.txt
else
tail -n +9 "$type_sugar"_"$sugar_ring"_"$ru"_scan.gjf > atoms_symbols.txt
fi

# Delete last two lines
head -n -2 atoms_symbols.txt > tmp.txt && mv tmp.txt atoms_symbols.txt 

# Print just first column of this file
awk '{ print $1 }' atoms_symbols.txt > tmp.txt && mv tmp.txt atoms_symbols.txt

# Print columns 4, 5, 6 of this file
awk '{ print  $1, "                 " }' atoms_symbols.txt > tmp.txt && mv tmp.txt atoms_symbols.txt

# Put all together 
paste atoms_symbols.txt cartesian_coord_min.txt > nucleoside.gjf

# Add charge and multiplicity
sed -i '1i0 1' nucleoside.gjf

# Add instructions
#sed -i '1 i\-1 1' nucleoside.gjf
echo " " > title_gjf_file.txt
echo "Opt and freq optimization for PES minimum (min#"$minimum_structure_number") of "$type_sugar"_"$sugar_ring"_"$ru" in "$environment" environment" >> title_gjf_file.txt
echo " " >> title_gjf_file.txt

if [ "$sugar_ring" == "none" ]; then

cat $MD_RUN/instructions_DFT/instructions_b3lyp.txt title_gjf_file.txt nucleoside.gjf > tmp.txt && mv tmp.txt "$type_sugar"_"$ru"_min.gjf
echo " " >> "$type_sugar"_"$ru"_min.gjf
# Replace name of chk file in final gjf file
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$type_sugar'_'$ru'_min.chk/g' "$type_sugar"_"$ru"_min.gjf
# Delete last 2 lines of this file
head -n -2 "$type_sugar"_"$ru"_min.gjf > tmp.txt && mv tmp.txt "$type_sugar"_"$ru"_min.gjf
# Substitute keywords for the case of water opt
if [ "$environment" == "water" ]; then
sed -i 's/scf=maxcycle=1600/scf=maxcycle=1600 scrf=(solvent=water)/g' "$type_sugar"_"$ru"_min.gjf
fi
# Replace instructions for DFT calculation
sed -i 's/opt=modredundant/opt=calcall freq/g' "$type_sugar"_"$ru"_min.gjf
sed -i 's/6-31g(d,p)/6-311++g(d,p)/g' "$type_sugar"_"$ru"_min.gjf
# Delete first line from file
sed -i 1,1d "$type_sugar"_"$ru"_min.gjf
# Replace RAM and PROC info
sed -i 's/%nprocshared=64/%nprocshared=16/g' "$type_sugar"_"$ru"_min.gjf
sed -i 's/%mem=160GB/%mem=16GB/g' "$type_sugar"_"$ru"_min.gjf
# Add empty line at the end
echo "" >> "$type_sugar"_"$ru"_min.gjf
# Copy file to corresponding folder
cp "$type_sugar"_"$ru"_min.gjf $MD_RUN/nucleosides/final_gaussian_opt/$environment/$type_sugar/$ru
# Copy job1 bash to this folder
cp $MD_RUN/job1.sh $MD_RUN/nucleosides/final_gaussian_opt/$environment/$type_sugar/$ru
# Replace name of file to calculate in job1
sed -i 's/g16 2deoxy_ribofuranose_2endo_alpha_adenine_scan.gjf/g16 '$type_sugar'_'$ru'_min.gjf/g' $MD_RUN/nucleosides/final_gaussian_opt/$environment/$type_sugar/$ru/job1.sh

else

cat $MD_RUN/instructions_DFT/instructions_b3lyp.txt title_gjf_file.txt nucleoside.gjf > tmp.txt && mv tmp.txt "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf
echo " " >> "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf
# Replace name of chk file in final gjf file
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$type_sugar'_'$sugar_ring'_'$ru'_min.chk/g' "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf
# Delete last 2 lines of this file
head -n -2 "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf > tmp.txt && mv tmp.txt "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf
# Substitute keywords for the case of water opt
if [ "$environment" == "water" ]; then
sed -i 's/scf=maxcycle=1600/scf=maxcycle=1600 scrf=(solvent=water)/g' "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf
fi
# Replace instructions for DFT calculation
sed -i 's/opt=modredundant/opt=calcall freq/g' "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf
sed -i 's/6-31g(d,p)/6-311++g(d,p)/g' "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf
# Delete first line from file
sed -i 1,1d "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf
# Replace RAM and PROC info
sed -i 's/%nprocshared=64/%nprocshared=16/g' "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf
sed -i 's/%mem=160GB/%mem=16GB/g' "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf
# Add empty line at the end
echo "" >> "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf
# Copy file to corresponding folder
cp "$type_sugar"_"$sugar_ring"_"$ru"_min.gjf $MD_RUN/nucleosides/final_gaussian_opt/$environment/$type_sugar/$sugar_ring/$ru
# Copy job1 bash to this folder
cp $MD_RUN/job1.sh $MD_RUN/nucleosides/final_gaussian_opt/$environment/$type_sugar/$sugar_ring/$ru
# Replace name of file to calculate in job1
sed -i 's/g16 2deoxy_ribofuranose_2endo_alpha_adenine_scan.gjf/g16 '$type_sugar'_'$sugar_ring'_'$ru'_min.gjf/g' $MD_RUN/nucleosides/final_gaussian_opt/$environment/$type_sugar/$sugar_ring/$ru/job1.sh

fi


rm *.txt


# THIS IS THE CASE THAT THE CALCULATION DID NOT FINISHED NORMALLY
else 






# Find if there is a line that says calculation could not be achived
if [ "$sugar_ring" == "none" ]; then
grep 'This type of calculation cannot be archived.' "$type_sugar"_"$ru"_scan2.log > message_about_archived_data.txt
else
grep 'This type of calculation cannot be archived.' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.log > message_about_archived_data.txt
fi
 
message_archived=`cat message_about_archived_data.txt` 
 
if [ "$message_archived" == " This type of calculation cannot be archived." ]; then


# Print angles values in different file
#for angle_value in $(seq 0 60 360); do
#echo $angle_value >> list_torsion_values.dat
#done

# Obtain list of relative energies
for structure_from_scan in $(seq 1 1 $num_conf_converged); do

# Get line number of each instance of -- Stationary point found.
if [ "$sugar_ring" == "none" ]; then
awk '/-- Stationary point found./{c++} c=='$structure_from_scan'{print NR;exit}' "$type_sugar"_"$ru"_scan2.log > line_number_conv_message.txt
line_number_with_conv_message=`cat line_number_conv_message.txt`
else
awk '/-- Stationary point found./{c++} c=='$structure_from_scan'{print NR;exit}' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.log > line_number_conv_message.txt
line_number_with_conv_message=`cat line_number_conv_message.txt`
fi

# Print line number for next convergence message -- Stationary point found.
if [ "$structure_from_scan" -eq 1 ]; then

# Copy from beginning to first convergence message -- Stationary point found.
if [ "$sugar_ring" == "none" ]; then
sed -n '1,'$line_number_with_conv_message'p' "$type_sugar"_"$ru"_scan2.log > structure#"$structure_from_scan"_file.txt 
else
sed -n '1,'$line_number_with_conv_message'p' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.log > structure#"$structure_from_scan"_file.txt 
fi

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
if [ "$sugar_ring" == "none" ]; then
awk '/-- Stationary point found./{c++} c=='$previous_structure_from_scan'{print NR;exit}' "$type_sugar"_"$ru"_scan2.log > previous_line_number_conv_message.txt
previous_line_number_with_conv_message=`cat previous_line_number_conv_message.txt`
sed -n ''$previous_line_number_with_conv_message','$line_number_with_conv_message'p' "$type_sugar"_"$ru"_scan2.log > structure#"$structure_from_scan"_file.txt 
else
awk '/-- Stationary point found./{c++} c=='$previous_structure_from_scan'{print NR;exit}' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.log > previous_line_number_conv_message.txt
previous_line_number_with_conv_message=`cat previous_line_number_conv_message.txt`
sed -n ''$previous_line_number_with_conv_message','$line_number_with_conv_message'p' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.log > structure#"$structure_from_scan"_file.txt 
fi

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

rm *.txt


# Print angles values in different file
#for angle_value in $(seq 0 60 360); do
#echo $angle_value >> list_torsion_values.dat
#done
 
 
fi





fi





