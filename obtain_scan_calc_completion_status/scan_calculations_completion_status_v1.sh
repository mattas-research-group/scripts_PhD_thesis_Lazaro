#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# Create folder for gaussian opt + freq of nucleosides
mkdir $MD_RUN/nucleosides/final_gaussian_opt
mkdir $MD_RUN/nucleosides/final_gaussian_opt/vacuum
mkdir $MD_RUN/nucleosides/final_gaussian_opt/water

function completion_status {

# Find all instances of convergence message and print how many are. If they are 7 then continue otherwise generate a new file to rerun the calculation
if [ "$sugar_ring" == "none" ]; then

grep -o 'Stationary point found' "$type_sugar"_"$ru"_scan.log | wc -l > number_conformers_that_converged.txt

else

grep -o 'Stationary point found' "$type_sugar"_"$sugar_ring"_"$ru"_scan.log | wc -l > number_conformers_that_converged.txt

fi

# Assign variable to this file
num_conf_converged=`cat number_conformers_that_converged.txt`

if [ "$num_conf_converged" == 7 ]; then

echo "The gaussian scan of the nucleoside of $type_sugar with $sugar_ring sugar ring conformation and RU $ru in $environment environment finished normally" >> $MD_RUN/list_calc_that_finished_normally.log

# Print angles values in different file
for angle_value in $(seq 0 60 360); do
echo $angle_value >> list_torsion_values.dat
done

# Now work on obtain the energies for each scan stationary point and their respective geometries and get the minimum
for structure_from_scan in $(seq 1 1 7); do

# Get line number of each instance of -- Stationary point found.
if [ "$sugar_ring" == "none" ]; then
awk '/-- Stationary point found./{c++} c=='$structure_from_scan'{print NR;exit}' "$type_sugar"_"$ru"_scan.log > line_number_conv_message.txt
line_number_with_conv_message=`cat line_number_conv_message.txt`
else
awk '/-- Stationary point found./{c++} c=='$structure_from_scan'{print NR;exit}' "$type_sugar"_"$sugar_ring"_"$ru"_scan.log > line_number_conv_message.txt
line_number_with_conv_message=`cat line_number_conv_message.txt`
fi

# Print line number for next convergence message -- Stationary point found.
if [ "$structure_from_scan" -eq 1 ]; then

# Copy from beginning to first convergence message -- Stationary point found.
if [ "$sugar_ring" == "none" ]; then
sed -n '1,'$line_number_with_conv_message'p' "$type_sugar"_"$ru"_scan.log > structure#"$structure_from_scan"_file.txt 
else
sed -n '1,'$line_number_with_conv_message'p' "$type_sugar"_"$sugar_ring"_"$ru"_scan.log > structure#"$structure_from_scan"_file.txt 
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
awk '/-- Stationary point found./{c++} c=='$previous_structure_from_scan'{print NR;exit}' "$type_sugar"_"$ru"_scan.log > previous_line_number_conv_message.txt
previous_line_number_with_conv_message=`cat previous_line_number_conv_message.txt`
sed -n ''$previous_line_number_with_conv_message','$line_number_with_conv_message'p' "$type_sugar"_"$ru"_scan.log > structure#"$structure_from_scan"_file.txt 
else
awk '/-- Stationary point found./{c++} c=='$previous_structure_from_scan'{print NR;exit}' "$type_sugar"_"$sugar_ring"_"$ru"_scan.log > previous_line_number_conv_message.txt
previous_line_number_with_conv_message=`cat previous_line_number_conv_message.txt`
sed -n ''$previous_line_number_with_conv_message','$line_number_with_conv_message'p' "$type_sugar"_"$sugar_ring"_"$ru"_scan.log > structure#"$structure_from_scan"_file.txt 
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
tail -n +9 *.gjf > atoms_symbols.txt

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



else # THIS IS THE CASE THAT THE CALCULATION DID NOT FINISHED NORMALLY






# Find if there is a line that says calculation could not be achived
if [ "$sugar_ring" == "none" ]; then
grep 'This type of calculation cannot be archived.' "$type_sugar"_"$ru"_scan.log > message_about_archived_data.txt
else
grep 'This type of calculation cannot be archived.' "$type_sugar"_"$sugar_ring"_"$ru"_scan.log > message_about_archived_data.txt
fi
 
message_archived=`cat message_about_archived_data.txt` 
 
if [ "$message_archived" == " This type of calculation cannot be archived." ]; then

echo "The gaussian scan of $type_sugar with $sugar_ring sugar ring conformation and RU $ru in $environment environment finished normally (the type of calculation cannot be archived)" >> $MD_RUN/list_calc_no_archived.log

# Print angles values in different file
for angle_value in $(seq 0 60 360); do
echo $angle_value >> list_torsion_values.dat
done

# Obtain list of relative energies
for structure_from_scan in $(seq 1 1 $num_conf_converged); do

# Get line number of each instance of -- Stationary point found.
if [ "$sugar_ring" == "none" ]; then
awk '/-- Stationary point found./{c++} c=='$structure_from_scan'{print NR;exit}' "$type_sugar"_"$ru"_scan.log > line_number_conv_message.txt
line_number_with_conv_message=`cat line_number_conv_message.txt`
else
awk '/-- Stationary point found./{c++} c=='$structure_from_scan'{print NR;exit}' "$type_sugar"_"$sugar_ring"_"$ru"_scan.log > line_number_conv_message.txt
line_number_with_conv_message=`cat line_number_conv_message.txt`
fi

# Print line number for next convergence message -- Stationary point found.
if [ "$structure_from_scan" -eq 1 ]; then

# Copy from beginning to first convergence message -- Stationary point found.
if [ "$sugar_ring" == "none" ]; then
sed -n '1,'$line_number_with_conv_message'p' "$type_sugar"_"$ru"_scan.log > structure#"$structure_from_scan"_file.txt 
else
sed -n '1,'$line_number_with_conv_message'p' "$type_sugar"_"$sugar_ring"_"$ru"_scan.log > structure#"$structure_from_scan"_file.txt 
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
awk '/-- Stationary point found./{c++} c=='$previous_structure_from_scan'{print NR;exit}' "$type_sugar"_"$ru"_scan.log > previous_line_number_conv_message.txt
previous_line_number_with_conv_message=`cat previous_line_number_conv_message.txt`
sed -n ''$previous_line_number_with_conv_message','$line_number_with_conv_message'p' "$type_sugar"_"$ru"_scan.log > structure#"$structure_from_scan"_file.txt 
else
awk '/-- Stationary point found./{c++} c=='$previous_structure_from_scan'{print NR;exit}' "$type_sugar"_"$sugar_ring"_"$ru"_scan.log > previous_line_number_conv_message.txt
previous_line_number_with_conv_message=`cat previous_line_number_conv_message.txt`
sed -n ''$previous_line_number_with_conv_message','$line_number_with_conv_message'p' "$type_sugar"_"$sugar_ring"_"$ru"_scan.log > structure#"$structure_from_scan"_file.txt 
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



else

 
 
echo "The gaussian scan of $type_sugar with $sugar_ring sugar ring conformation and RU $ru in $environment environment DID NOT FINISH NORMALLY" >> $MD_RUN/list_calc_that_didnt_finished_normally.log 

# Print angles values in different file
for angle_value in $(seq 0 60 360); do
echo $angle_value >> list_torsion_values.dat
done
 
 
 
fi

# Generate second set of calculations to finish gaussian scan normally
# Convert log file to gjf file
if [ "$sugar_ring" == "none" ]; then
obabel "$type_sugar"_"$ru"_scan.log -O nucleoside.gjf
else
obabel "$type_sugar"_"$sugar_ring"_"$ru"_scan.log -O nucleoside.gjf
fi

# Delete first two lines from final gjf file
sed -i 1,2d nucleoside.gjf

# Attach DFT instructions
cat $MD_RUN/instructions_DFT/instructions_b3lyp.txt nucleoside.gjf > nucleoside2.gjf

# Delete instructions for cherif cluster
sed -i 1,1d nucleoside2.gjf

# Substitute keywords for the case of water opt
if [ "$environment" == "water" ]; then
sed -i 's/scf=maxcycle=1600/scf=maxcycle=1600 scrf=(solvent=water)/g' nucleoside2.gjf
fi

# Change name of file and attach job1.sh with correct name
if [ "$sugar_ring" == "none" ]; then

# Change name of file
mv nucleoside2.gjf "$type_sugar"_"$ru"_scan2.gjf

# Copy in current folder job1
cp $MD_RUN/job1.sh .

# Replace name of file by backbone name in job1 script
sed -i 's/g16 2deoxy_ribofuranose_2endo_alpha_adenine_scan.gjf/g16 '$type_sugar'_'$ru'_scan2.gjf/g' job1.sh

# Replace name of chk file in second scan file
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$type_sugar'_'$ru'_scan.chk/g' "$type_sugar"_"$ru"_scan2.gjf

# Replace gaussian instructions for recalculate the scan
sed -i 's/opt=modredundant/opt=(tight,readfc) geom=allcheck guess=read int=ultrafine/g' "$type_sugar"_"$ru"_scan2.gjf

# Replace RAM and PROC info
sed -i 's/%nprocshared=64/%nprocshared=16/g' "$type_sugar"_"$ru"_scan2.gjf
sed -i 's/%mem=160GB/%mem=16GB/g' "$type_sugar"_"$ru"_scan2.gjf

else

# Rename file
mv nucleoside2.gjf "$type_sugar"_"$sugar_ring"_"$ru"_scan2.gjf

# Copy job script in current folder
cp $MD_RUN/job1.sh .

# Replace name of file by backbone name
sed -i 's/g16 2deoxy_ribofuranose_2endo_alpha_adenine_scan.gjf/g16 '$type_sugar'_'$sugar_ring'_'$ru'_scan2.gjf/g' job1.sh

# Replace name of chk file in second scan file
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$type_sugar'_'$sugar_ring'_'$ru'_scan.chk/g' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.gjf

# Replace gaussian instructions for recalculate the scan
sed -i 's/opt=modredundant/opt=(tight,readfc) geom=allcheck guess=read int=ultrafine/g' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.gjf

# Replace RAM and PROC info
sed -i 's/%nprocshared=64/%nprocshared=16/g' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.gjf
sed -i 's/%mem=160GB/%mem=16GB/g' "$type_sugar"_"$sugar_ring"_"$ru"_scan2.gjf

fi



fi

}

##########################################################################################################################
##########################################################################################################################

############################################## WORK ON EACH FOLDER #######################################################

# Open corresponding folder
cd nucleosides/gaussian_scan

echo "vacuum" > environment.txt
echo "water" >> environment.txt

cat environment.txt | while read environment

do

cd $environment

cat $MD_RUN/list_folders/sugars_list.txt | while read sugars

do

if [[ "$sugars" == *\/* ]]; then # This is the case that the TC is a sugar with sugar ring conformation

# Obtain trifunctional connector name
echo $sugars | cut -d "/" -f1 > sugar_name.txt
type_sugar=`cat sugar_name.txt`

# Obtain sugar ring conformation name
echo $sugars | cut -d "/" -f2 > sugar_ring_conf.txt
sugar_ring=`cat sugar_ring_conf.txt`

# Create corresponding folder in final gaussian folder
mkdir $MD_RUN/nucleosides/final_gaussian_opt/$environment/$type_sugar
mkdir $MD_RUN/nucleosides/final_gaussian_opt/$environment/$type_sugar/$sugar_ring

cd $type_sugar/$sugar_ring

cat $MD_RUN/list_folders/bases_list.txt | while read ru

do

# Create corresponding folder in final gaussian folder
mkdir $MD_RUN/nucleosides/final_gaussian_opt/$environment/$type_sugar/$sugar_ring/$ru

cd $ru

# Message
echo "working on $environment environment, $type_sugar with $sugar_ring sugar ring conformation and base $ru"

# Create file that will contain list of angles, energies and sugar ring conf 
touch sugar_ring_conformation_name_repeated.txt
touch list_torsion_values.txt
#touch list_relative_energies_kcal_mol.txt
touch RU_name_repeated.txt

# Print outside file with list of sugar ring conformations
for i in {1..7}; do
echo "$sugar_ring" >> sugar_ring_conformation_name_repeated.txt
done 

# Print outside file with list of sugar ring conformations
for i in {1..7}; do
echo "$ru" >> RU_name_repeated.txt
done

# Run function that determines if calculation completed
completion_status













cd ../

done

cd ../../

else # This is the case that the TC is glycerol or glyceric acid or peptide which is not a sugar and does not has sugar ring conformation

# Obtain trifunctional connector name
echo $sugars > sugar_name.txt
type_sugar=`cat sugar_name.txt`

# Obtain sugar ring conformation name
sugar_ring="none"

# Create corresponding folder in final gaussian folder
mkdir $MD_RUN/nucleosides/final_gaussian_opt/$environment/$type_sugar

if [ "$type_sugar" != "peptide" ]; then

cd $type_sugar

cat $MD_RUN/list_folders/bases_list.txt | while read ru

do

# Create corresponding folder in final gaussian folder
mkdir $MD_RUN/nucleosides/final_gaussian_opt/$environment/$type_sugar/$ru

cd $ru

sugar_ring="none"

# Message
echo "working on $environment environment, $type_sugar with $sugar_ring sugar ring conformation and base $ru"

# Create file that will contain list of angles, energies and sugar ring conf 
touch list_torsion_values.txt
#touch list_relative_energies_kcal_mol.txt
touch RU_name_repeated.txt

# Print outside file with list of sugar ring conformations
for i in {1..7}; do
echo "$ru" >> RU_name_repeated.txt
done

# Run function that determines if calculation completed
completion_status










cd ../

done

elif [ "$type_sugar" == "peptide" ]; then

cd $type_sugar

cat $MD_RUN/list_folders/bases_ac_list.txt | while read ru

do

# Create corresponding folder in final gaussian folder
mkdir $MD_RUN/nucleosides/final_gaussian_opt/$environment/$type_sugar/$ru

cd $ru

sugar_ring="none"

# Run function that determines if calculation completed
completion_status













cd ../

done

fi


cd ../

fi

done

cd ../

done

