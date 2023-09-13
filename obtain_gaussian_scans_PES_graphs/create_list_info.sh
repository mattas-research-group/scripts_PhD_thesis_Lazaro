#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)


############################################## WORK ON EACH FOLDER #######################################################

# Open corresponding folder
cd nucleosides/gaussian_scan

declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid" "peptide")

declare -a environment=("vacuum" "water")

for environment in "${environment[@]}"; do

cd $environment

for TC in "${TC_names[@]}"; do

cd $TC

# Create general files that contain list of names for ring conformation
touch list_Ring_Conf_names.txt
echo "Ring_Conformation" > list_Ring_Conf_names.txt

# Create general files that contain list of names for RU
touch list_RU_names.txt
echo "Recognition_Unit" > list_RU_names.txt

# Create general files that contain list of relative energies
touch list_rEnergies.txt
echo "Relative_Energies(Kcal/mol)" > list_rEnergies.txt

# Create general files that contain list of torsion angles values in degrees
touch list_torsion_angles.txt
echo "Torsion_Angle(degrees)" > list_torsion_angles.txt

cd ../
done
cd ../
done





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





cd $type_sugar/$sugar_ring

cat $MD_RUN/list_folders/bases_list.txt | while read ru

do

cd $ru

# Message
#echo "working on $environment environment, $type_sugar with $sugar_ring sugar ring conformation and base $ru"

# Print information in general files
for s in {1..7}; do
echo "$sugar_ring"
echo "$sugar_ring" >> $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/list_Ring_Conf_names.txt
done

for r in {1..7}; do
echo "$ru" >> $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/list_RU_names.txt
done

list_energies=`cat list_relative_energies_scan_kcal_mol.dat`
echo "$list_energies" >> $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/list_rEnergies.txt

list_torsions=`cat list_torsion_values.dat`
echo "$list_torsions" >> $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/list_torsion_angles.txt


cd ../

done

cd ../../


else # This is the case that the TC is glycerol or glyceric acid or peptide which is not a sugar and does not has sugar ring conformation

# Obtain trifunctional connector name
echo $sugars > sugar_name.txt
type_sugar=`cat sugar_name.txt`

# Obtain sugar ring conformation name
sugar_ring="none"



if [ "$type_sugar" != "peptide" ]; then

cd $type_sugar

cat $MD_RUN/list_folders/bases_list.txt | while read ru

do

cd $ru

sugar_ring="none"

# Print information in general files
for r in {1..7}; do
echo "$sugar_ring" >> $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/list_Ring_Conf_names.txt
done

for r in {1..7}; do
echo "$ru" >> $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/list_RU_names.txt
done

list_energies=`cat list_relative_energies_scan_kcal_mol.dat`
echo "$list_energies" >> $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/list_rEnergies.txt

list_torsions=`cat list_torsion_values.dat`
echo "$list_torsions" >> $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/list_torsion_angles.txt






cd ../

done

elif [ "$type_sugar" == "peptide" ]; then

cd $type_sugar

cat $MD_RUN/list_folders/bases_ac_list.txt | while read ru

do

cd $ru

sugar_ring="none"

# Print information in general files
for r in {1..7}; do
echo "$sugar_ring" >> $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/list_Ring_Conf_names.txt
done

for r in {1..7}; do
echo "$ru" >> $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/list_RU_names.txt
done

list_energies=`cat list_relative_energies_scan_kcal_mol.dat`
echo "$list_energies" >> $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/list_rEnergies.txt

list_torsions=`cat list_torsion_values.dat`
echo "$list_torsions" >> $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/list_torsion_angles.txt







cd ../

done

fi


cd ../

fi

done

cd ../

done
