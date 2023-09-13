#!/bin/bash

# Note1: no_int_coord_base.gzmat and no_int_coord_sugar.gzmat are the files with only the z matrix
# Note2: this initial test was with beta 2deoxyribofuranose and guanine

# Insert the names of the files for the TC and the RUs
read -p "insert name of log file for the Trifunctional Connector (TC) without extension:" tc
read -p "insert name of log file for the Recognition unit (RU) without extension:" ru

# Convert both files for the sugar and base to gzmat
obabel *.log -O *.gzmat

########################
## WORK ON SUGAR FILE ##
########################

# Declare labels of anomeric carbon
read -p "insert label of anomeric carbon as appears in gaussview:" anomer ####### ************************(THIS LINE CAN CHANGE DEPENDING ON THE TC)***********************

# Declare label of hemiacetalic oxygen for beta or 
read -p "insert label of hemiacetalic oxygen as appears in gaussview:" hemiacetal_O ####### ************************(THIS LINE CAN CHANGE DEPENDING ON THE TC)***********************

# First declare the labels of the atoms to be deleted from the TC (in sugar case should be an oxygen and hydrogen)
read -p "insert labels of oxygen and hydrogen (OH) to be deleted as appears in gaussview(e.g. 17 18):" oxygen hydrogen ########## ************************(THIS LINE CAN CHANGE DEPENDING ON THE TC)***********************

# Obtain the lines number for these atoms
line_oxygen=$(($oxygen+6))
line_hydrogen=$(($hydrogen+6))

# Delete lines in matrix that correspond to the O and H in anomeric position
sed ''$line_oxygen','$line_hydrogen'd' $tc.gzmat > sugar_sin_OH.gzmat

#### Delete lines with values of the distances, angles and dihedrals for the atoms deleted in previous step ###
sed -i '/r'$oxygen'/d' sugar_sin_OH.gzmat
sed -i '/a'$oxygen'/d' sugar_sin_OH.gzmat
sed -i '/d'$oxygen'/d' sugar_sin_OH.gzmat
sed -i '/r'$hydrogen'/d' sugar_sin_OH.gzmat
sed -i '/a'$hydrogen'/d' sugar_sin_OH.gzmat
sed -i '/d'$hydrogen'/d' sugar_sin_OH.gzmat

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
lines=$(($atoms_sugar_ion+7))

################################################################

# Substitute the last r# a# and d# by number of atoms in sugar ion
sed -i -e 's/r'$atoms_sugar'/r'$atoms_sugar_ion'/g' sugar_sin_OH.gzmat
sed -i -e 's/a'$atoms_sugar'/a'$atoms_sugar_ion'/g' sugar_sin_OH.gzmat
sed -i -e 's/d'$atoms_sugar'/d'$atoms_sugar_ion'/g' sugar_sin_OH.gzmat

# Copy lines with values of internal coordinates to another file
sed -n ''$lines',$p' sugar_sin_OH.gzmat > internal_coordinates_values_sugar_sinOH.txt 
# Eliminate last empty line
sed -i '$ d' internal_coordinates_values_sugar_sinOH.txt      ### THIS IS THE FINAL FILE FOR INTERNAL COORDINATES VALUES OF THE SUGAR

# Delete values the internal coordinates
sed -i '1,/Variables:/!d' sugar_sin_OH.gzmat

# Delete last line of this file 
sed -i '$d' sugar_sin_OH.gzmat

# Change the name
mv sugar_sin_OH.gzmat no_int_coord_sugar_sinOH.gzmat ### THIS IS THE FINAL FILE FOR THE Z MATRIX OF SUGAR WITHOUT THE INT COORD VALUES

rm sinheadings_sugar.gzmat number_atoms_sugar.txt 


########################
## WORK ON BASE FILE  ##
########################

# Delete headings from base file
sed -i '1,6d' $ru.gzmat

### Calculate number of atoms in base #######

# Delete values of the internal coordinates
sed '1,/Variables:/!d' $ru.gzmat > no_int_coord_base.txt

# Delete last line of this file 
sed -i '$d' no_int_coord_base.txt

# Count number of lines in this file and determine number of atoms in base
< no_int_coord_base.txt wc -l > number_atoms_base_new.txt

atoms_base=`cat number_atoms_base_new.txt`
echo $atoms_base
ions_atoms_base=$(($atoms_base-1))
echo $ions_atoms_base
first_atom_base=$(($atoms_sugar_ion+1))
second_atom_base=$(($atoms_sugar_ion+2))
third_atom_base=$(($atoms_sugar_ion+3))

# Delete last line to create z matrix for ionic base: same as delete hydrogen from base to create ion 
read -p "insert label of hydrogen atom to be deleted from base as appears in gaussview:" hydrogen ####### ************************(THIS LINE CAN CHANGE DEPENDING ON THE RU)***********************
sed -i ''$hydrogen'd' no_int_coord_base.txt

# Change labels for distance, angle and dihedrals in this file
labels_number_base=$(($atoms_sugar_ion+4))
final_labels_base=$(($atoms_sugar_ion+$ions_atoms_base))
echo $labels_number_base
echo $final_labels_base 

### Change r#, a# and d# and atoms labels in file of base from line 4 by new numbers

for lines in $(seq 4 1 $ions_atoms_base)

do

# Calculate new r# a# and d#
new_labels_number_base=$(($atoms_sugar_ion+$lines))
sed -i ''$lines's/r'$lines'/r'$new_labels_number_base'/' no_int_coord_base.txt
sed -i ''$lines's/a'$lines'/a'$new_labels_number_base'/' no_int_coord_base.txt
sed -i ''$lines's/d'$lines'/d'$new_labels_number_base'/' no_int_coord_base.txt

# Calculate new atoms labels for the 
awk 'NR=='$lines' { print $1, $2=$2+'$atoms_sugar_ion', $3, $4=$4+'$atoms_sugar_ion', $5, $6=$6+'$atoms_sugar_ion', $7 }' no_int_coord_base.txt >> no_int_coord_base_labels_fixed.txt 

done

## Add the three first lines to the z matriz of the base and fix them accordingly

# Replace first line with internal coordinates as it should be in file
sed -i '1 s/.*/\N  '$anomer'  r'$first_atom_base'  '$hemiacetal_O'  a'$first_atom_base'  3  d'$first_atom_base'/' no_int_coord_base.txt

# Replace second line with another specific coordinates
sed -i '2 s/.*/\C  '$first_atom_base'  r'$second_atom_base'  '$anomer'  a'$second_atom_base'  '$hemiacetal_O'  d'$second_atom_base'/' no_int_coord_base.txt

# Replace third line with another specific coordinates
sed -i '3 s/.*/\C  '$second_atom_base'  r'$third_atom_base'  '$first_atom_base'  a'$third_atom_base'  '$anomer'  d'$third_atom_base'/' no_int_coord_base.txt

# Delete lines from line 3 in no_int_coord_base.txt
sed -i '4,$d' no_int_coord_base.txt

# Paste all together 
cat no_int_coord_base.txt no_int_coord_base_labels_fixed.txt > no_int_coord_base_labels_fixed_v2.txt

######## Replace first hydrogen line with correct atoms labels #######

# First output number of line with first H
awk '/H/{ print NR; exit }' no_int_coord_base_labels_fixed_v2.txt > line_of_first_H.txt
line_H=`cat line_of_first_H.txt`
echo $line_H

for lines in $(seq 1 1 $ions_atoms_base)

do

if [[ $lines -eq $line_H ]]
then
# Change the actual line
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else 
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi

done

rm no_int_coord_base_labels_fixed_v2.txt no_int_coord_base_labels_fixed.txt no_int_coord_base.txt

# Rename this file as final
mv no_int_coord_base_labels_fixed_v3.txt no_int_coord_base_sinH.gzmat  #### THIS IS THE FINAL FILE FOR THE Z MATRIZ OF THE BASE FIXED APPROPIATELY



################ Work in internal coordinates values of base ############################

# Copy values of internal coordinates in different file
lines2=$(($ions_atoms_base+1))

# Copy lines with values of internal coordinates to another file
sed -n '/Variables:/,$p' $ru.gzmat > internal_coordinates_values_base.txt

# Delete first line of this file
sed -i '1d' internal_coordinates_values_base.txt

# Add new lines with new values of distances, angles and dihedrals to file with internal coordinates (this is where I add the numbers of internal coordinates to keep in all molecules)*****IMPORTANT!

echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 109.59" >> first_int_coord_base.txt
echo "d1= 232.30" >> first_int_coord_base.txt

# Insert a2 d2 and d3 in this file 
sed -i '2 i a2= 131.20' internal_coordinates_values_base.txt
sed -i '3 i d2= 81.59' internal_coordinates_values_base.txt
sed -i '6 i d3= 177.00' internal_coordinates_values_base.txt

# Put both files together
cat first_int_coord_base.txt internal_coordinates_values_base.txt > internal_coordinates_base.txt

# Delete coordinates values for hydrogen that was deleted from base

hydrogen_deleted=16 # This have to be modify accordingly ************************(THIS LINE CAN CHANGE DEPENDING ON THE RU)***********************

sed -i '/r'$hydrogen_deleted'/d' internal_coordinates_base.txt
sed -i '/a'$hydrogen_deleted'/d' internal_coordinates_base.txt
sed -i '/d'$hydrogen_deleted'/d' internal_coordinates_base.txt

# Change labels of each atoms in internal coordinates file of base

for labels in $(seq 1 1 $ions_atoms_base)

do

new_label=$(($labels+$atoms_sugar_ion))

echo "r$new_label=" >> newlabels_base.txt
echo "a$new_label=" >> newlabels_base.txt
echo "d$new_label=" >> newlabels_base.txt

done

# Eliminate everything before the = sign in old int coord file from the base
cut -f2 -d"=" internal_coordinates_base.txt > new_internal_coordinates_base.txt

# Paste new labels and values for internal coordinates
paste newlabels_base.txt new_internal_coordinates_base.txt > new_internal_coordinates_base2.txt

# Reduce spaces between two columns
sed -i $'s/\t//g' new_internal_coordinates_base2.txt

tr -s " " < new_internal_coordinates_base2.txt > internal_coordinates_base_sinH.gzmat #### THIS IS THE FINAL FILE WITH THE INTERNAL COORD VALUES FIXED FOR THE BASE 






# Fix value of dihedral angle between first hydrogen atom of base 
# First print dihedral angle appropiately
awk '(NR=='$line_H'){print $7}' no_int_coord_base_sinH.gzmat > dihedral_level_firstH.txt
value=`cat dihedral_level_firstH.txt`

# Output number of line in internal coordinate file with this value
grep -n "$value" internal_coordinates_base_sinH.gzmat | cut -d : -f1 > line_number_in_intcoord_base.txt

line=`cat line_number_in_intcoord_base.txt`

# Replace the line
awk 'NR=='$line' {$0="'$value'= 4.26"} 1' internal_coordinates_base_sinH.gzmat > internal_coordinates_base_sinH_v2.gzmat






################### Now paste all together to create final z matriz file for the nucleosides ######################

cat no_int_coord_sugar_sinOH.gzmat no_int_coord_base_sinH.gzmat internal_coordinates_values_sugar_sinOH.txt internal_coordinates_base_sinH_v2.gzmat > nucleoside.gzmat

# Add instructions for calculations
read -p "insert name of your output file gjf:" nucleoside

# Convert this file to gjf
obabel nucleoside.gzmat -O $nucleoside.gjf

# Paste instructions for b3lyp
sed -i -e 1,2d $nucleoside.gjf
cat instructions_DFT/instructions_b3lyp.txt $nucleoside.gjf > "$nucleoside"_new.gjf

# Replace name chk by name of cell
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$nucleoside'.chk/g' "$nucleoside"_new.gjf

# Add frozen coordinates
cat "$nucleoside"_new.gjf frozen_coordinates/frozen_coord_fura.txt > "$nucleoside"_new2.gjf

rm $nucleoside.gjf "$nucleoside"_new.gjf

mv "$nucleoside"_new2.gjf $nucleoside.gjf



