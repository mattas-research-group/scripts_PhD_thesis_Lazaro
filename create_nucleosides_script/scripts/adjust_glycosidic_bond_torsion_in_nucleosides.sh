#!/bin/bash

# Insert the names of the files for the TC and the RUs
read -p "insert name of log file for the Recognition unit (RU) without extension:" ru
#read -p "insert type of anomer in case that it is (e.g. beta, or alpha, or none if there is not anomer)" type_anomer
read -p "insert type of trifunctional connector (TC) in case that it is (e.g. 2deoxy_ribofuranose, or ribofuranose, or 2deoxy_ribopyranose, or ribopyranose, or threose, or glycerol, or peptide)" type_sugar
read -p "insert conformation of the sugar ring (e.g. 2endo_alpha, or 2endo_beta, or D1C4_alpha, or D4C1_beta, or none if it is glycerol, glyceric acid or peptide)" sugar_ring

################ Insert options to authomatically change labels of atoms to delete and hemiacetal O and anomeric carbon for sugar ########################

if [ "$type_sugar" == "2deoxy_ribofuranose" ]; then

anomer=8
hemiacetal_O=4

elif [ "$type_sugar" == "ribofuranose" ]; then

anomer=8
hemiacetal_O=4

elif [ "$type_sugar" == "threose" ]; then

anomer=5
hemiacetal_O=2

elif [ "$type_sugar" == "2deoxy_ribopyranose" ]; then

anomer=1
hemiacetal_O=2

elif [ "$type_sugar" == "ribopyranose" ]; then

anomer=1
hemiacetal_O=2

elif [ "$type_sugar" == "glycerol" ]; then

anomer=1
hemiacetal_O=2

elif [ "$type_sugar" == "glyceric_acid" ]; then
anomer=10
hemiacetal_O=1

elif [ "$type_sugar" == "peptide" ]; then
anomer=1
hemiacetal_O=2

fi



# Give name to trifunctional connector file
if [ "$type_sugar" == "glycerol" ] || [ "$type_sugar" == "glyceric_acid" ] || [ "$type_sugar" == "peptide" ]; then
tc=$type_sugar
else
tc="$type_sugar"_"$sugar_ring"
fi

# Convert nucleoside log file to gzmat
nucleoside="$tc"_"$ru"
obabel $nucleoside.log -O $nucleoside.gzmat

# Convert first TC log or mol file to gzmat
# Convert tc to gzmat
if [ "$sugar_ring" == "D1C4_beta" ] || [ "$sugar_ring" == "D4C1_beta" ] || [ "$type_sugar" == "peptide" ]; then
obabel $tc.mol -O $tc.gzmat
else
obabel $tc.log -O $tc.gzmat
fi

# Eliminate lines to leave only z matrix minus two atoms for the TC
# Delete first headings of this file
sed '1,6d' $tc.gzmat > sinheadings_sugar.gzmat

# Delete values the internal coordinates
sed -i '1,/Variables:/!d' sinheadings_sugar.gzmat
sed -i '$d' sinheadings_sugar.gzmat

# Count number of lines in TC gzmat file without headings and int coord values 
< sinheadings_sugar.gzmat wc -l > number_atoms_sugar.txt
atoms_sugar=`cat number_atoms_sugar.txt`
atoms_sugar_ion=$(($atoms_sugar-2))

first_base_atom=$(($atoms_sugar-1))
torsion_to_modif=$(($atoms_sugar_ion+2))

echo "the torsion to change correspond to atom#$torsion_to_modif from the base"

# Change this atom in int coord values of z matrix of nucleoside
if [ "$sugar_ring" == "D1C4_beta" ] || [ "$sugar_ring" == "D4C1_beta" ] || [ "$sugar_ring" == "2endo_beta" ] || [ "$sugar_ring" == "3endo_beta" ]; then

sed -i 's/d'$torsion_to_modif'=.*/d'$torsion_to_modif'= 325/' $nucleoside.gzmat # The value may change

elif [ "$sugar_ring" == "D1C4_alpha" ] || [ "$sugar_ring" == "D4C1_alpha" ] || [ "$sugar_ring" == "2endo_alpha" ] || [ "$sugar_ring" == "3endo_alpha" ]; then

sed -i 's/d'$torsion_to_modif'=.*/d'$torsion_to_modif'= 325/' $nucleoside.gzmat # The value may change

else

sed -i 's/d'$torsion_to_modif'=.*/d'$torsion_to_modif'= 325/' $nucleoside.gzmat # The value may change

fi

# Convert nucleoside new gzmat file to gjf
obabel $nucleoside.gzmat -O $nucleoside.gjf

# Paste instructions for b3lyp
sed -i -e 1,2d $nucleoside.gjf

cat instructions_DFT/instructions_b3lyp.txt $nucleoside.gjf > "$nucleoside"_new.gjf

# Replace name chk by name of cell
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$nucleoside'.chk/g' "$nucleoside"_new.gjf

# Paste at the end coordinates to freeze
if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ]; then
cat "$nucleoside"_new.gjf frozen_coordinates/frozen_coord_fura.txt > "$nucleoside"_new2.gjf
elif [ "$type_sugar" == "threose" ]; then
cat "$nucleoside"_new.gjf frozen_coordinates/frozen_coord_threo.txt > "$nucleoside"_new2.gjf
elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ]; then
cat "$nucleoside"_new.gjf frozen_coordinates/frozen_coord_pyra.txt > "$nucleoside"_new2.gjf
else
mv "$nucleoside"_new.gjf "$nucleoside"_new2.gjf
fi

# Add coordinates to be scanned
echo " " > coordinates_to_scan.txt
echo "* $first_base_atom $anomer *" >> coordinates_to_scan.txt
echo "D $torsion_to_modif $first_base_atom $anomer $hemiacetal_O S 12 30.000000" >> coordinates_to_scan.txt

cat "$nucleoside"_new2.gjf coordinates_to_scan.txt > "$nucleoside"_scan.gjf











