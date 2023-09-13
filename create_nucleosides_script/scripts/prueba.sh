#!/bin/bash

atoms_sugar_ion=17
ions_atoms_base=15

for lines in $(seq 4 1 $ions_atoms_base)

do

# Calculate new r# a# and d#
new_labels_number_base=$(($atoms_sugar_ion+$lines))
sed -i ''$lines's/r'$lines'/r'$new_labels_number_base'/' no_int_coord_base.txt
sed -i ''$lines's/a'$lines'/a'$new_labels_number_base'/' no_int_coord_base.txt
sed -i ''$lines's/d'$lines'/d'$new_labels_number_base'/' no_int_coord_base.txt

# Calculate new atoms labels
awk 'NR=='$lines' { print $1,  $2=$2+'$atoms_sugar_ion',  $3,  $4=$4+'$atoms_sugar_ion',  $5,  $6=$6+'$atoms_sugar_ion',  $7 }' no_int_coord_base.txt >> fixed.txt 

done









