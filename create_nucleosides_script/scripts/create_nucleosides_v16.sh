#!/bin/bash

# Note1: no_int_coord_base.gzmat and no_int_coord_sugar.gzmat are the files with only the z matrix
# Note2: this initial test was with beta 2deoxyribofuranose and guanine

# Insert the names of the files for the TC and the RUs
#read -p "insert name of log file for the Trifunctional Connector (TC) without extension:" tc
read -p "insert name of log file for the Recognition unit (RU) without extension:" ru
#read -p "insert type of anomer in case that it is (e.g. beta, or alpha, or none if there is not anomer)" type_anomer
read -p "insert type of sugar in case that it is (e.g. 2deoxy_ribofuranose, or ribofuranose, or 2deoxy_ribopyranose, or ribopyranose, or threose, or none if the TC is not a sugar)" type_sugar
read -p "insert conformation of the sugar ring (e.g. 2endo_alpha, or 2endo_beta, or D1C4_alpha, or D4C1_beta)" sugar_ring

# Give name to trifunctional connector file
tc="$type_sugar"_"$sugar_ring"

# Obtain type of anomer
cut -d "_" -f2 <<< "$sugar_ring" > type_anomer.txt
type_anomer=`cat type_anomer.txt`

# Convert ru to gzmat
if [ "$ru" == "BA_cbond" ] || [ "$ru" == "BA_nbond" ] || [ "$ru" == "TARC_cbond" ] || [ "$ru" == "TARC_nbond" ] || [ "$ru" == "melamine" ]; then
obabel $ru.mol -O $ru.gzmat
else
obabel $ru.log -O $ru.gzmat
fi

# Convert tc to gzmat
if [ "$sugar_ring" == "D1C4_beta" ] || [ "$sugar_ring" == "D4C1_beta" ]; then
obabel $tc.mol -O $tc.gzmat
else
obabel $tc.log -O $tc.gzmat
fi

################ Insert options to authomatically change labels of atoms to delete and hemiacetal O and anomeric carbon for sugar ########################

if [ "$type_sugar" == "2deoxy_ribofuranose" ]; then

anomer=8
hemiacetal_O=4
oxygen=17
hydrogen=18

elif [ "$type_sugar" == "ribofuranose" ]; then

anomer=8
hemiacetal_O=4
oxygen=19
hydrogen=20

elif [ "$type_sugar" == "threose" ]; then

anomer=5
hemiacetal_O=2
oxygen=15
hydrogen=16

elif [ "$type_sugar" == "2deoxy_ribopyranose" ]; then

anomer=1
hemiacetal_O=2
oxygen=17
hydrogen=18

elif [ "$type_sugar" == "ribopyranose" ]; then

anomer=1
hemiacetal_O=2
oxygen=19
hydrogen=20

fi


########################
## WORK ON SUGAR FILE ##
########################

# Declare labels of anomeric carbon
#read -p "insert label of anomeric carbon as appears in gaussview:" anomer ####### ************************(THIS LINE CAN CHANGE DEPENDING ON THE TC)***********************

# Declare label of hemiacetalic oxygen for beta or 
#read -p "insert label of hemiacetalic oxygen as appears in gaussview:" hemiacetal_O ####### ************************(THIS LINE CAN CHANGE DEPENDING ON THE TC)***********************

if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ]; then
atom_before_hemiacetal=$(($hemiacetal_O-1))
atom_before_before_hemiacetal=$(($hemiacetal_O-2))
elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ]; then
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
lines_plus_atoms_in_sugar_ion=$(($atoms_sugar_ion+6))
starting_line_int_coord=$(($atoms_sugar_ion+7))

################################################################

# Substitute the last r# a# and d# by number of atoms in sugar ion
#sed -i -e 's/r'$atoms_sugar'/r'$atoms_sugar_ion'/g' sugar_sin_OH.gzmat
#sed -i -e 's/a'$atoms_sugar'/a'$atoms_sugar_ion'/g' sugar_sin_OH.gzmat
#sed -i -e 's/d'$atoms_sugar'/d'$atoms_sugar_ion'/g' sugar_sin_OH.gzmat

# Copy lines with values of internal coordinates to another file
sed -n ''$starting_line_int_coord',$p' sugar_sin_OH.gzmat > internal_coordinates_values_sugar_sinOH.txt 
# Eliminate last empty line
sed -i '$ d' internal_coordinates_values_sugar_sinOH.txt      

# Delete values the internal coordinates
sed -i '1,/Variables:/!d' sugar_sin_OH.gzmat

# Delete last line of this file 
sed -i '$d' sugar_sin_OH.gzmat

############# Work on correcting the rlabel, alabel and dlabel in the file without internal coordinates values ################
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


########################
## WORK ON BASE FILE  ##
########################

################ Insert options to authomatically change labels of atoms to delete for base ########################

if [ "$ru" == "adenine" ] || [ "$ru" == "thymine" ]; then

hydrogen=15

elif [ "$ru" == "guanine" ]; then

hydrogen=16

elif [ "$ru" == "cytosine" ]; then

hydrogen=13

elif [ "$ru" == "uracil" ] || [ "$ru" == "CA" ]; then

hydrogen=12

elif [ "$ru" == "melamine" ]; then

hydrogen=15

elif [ "$ru" == "BA_cbond" ]; then

hydrogen=13

elif [ "$ru" == "BA_nbond" ]; then

hydrogen=11

elif [ "$ru" == "TARC_cbond" ] || [ "$ru" == "TARC_nbond" ]; then

hydrogen=16

fi


# Delete headings from base file
sed -i '1,6d' $ru.gzmat

### Calculate number of atoms in base #######

# Delete values the internal coordinates
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
fourth_atom_base=$(($atoms_sugar_ion+4))
six_atom_base=$(($atoms_sugar_ion+6))
seven_atom_base=$(($atoms_sugar_ion+7))
nine_atom_base=$(($atoms_sugar_ion+9))

# Delete last line to create z matrix for ionic base: same as delete hydrogen from base to create ion 
#read -p "insert label of hydrogen atom to be deleted from base as appears in gaussview:" hydrogen ####### ************************(THIS LINE CAN CHANGE DEPENDING ON THE RU)***********************
sed -i ''$hydrogen'd' no_int_coord_base.txt

# Change labels for distance, angle and dihedrals in this file
labels_number_base=$(($atoms_sugar_ion+4))
final_labels_base=$(($atoms_sugar_ion+$ions_atoms_base))
echo $labels_number_base
echo $final_labels_base 

### Change r#, a# and d# and atoms labels in file of base from line 4 by new numbers

for lines in $(seq 4 1 $ions_atoms_base)

do

if [[ "$lines" -lt "$hydrogen" ]]; then
# Calculate new r# a# and d#
new_labels_number_base=$(($atoms_sugar_ion+$lines))
sed -i ''$lines's/r'$lines'/r'$new_labels_number_base'/' no_int_coord_base.txt
sed -i ''$lines's/a'$lines'/a'$new_labels_number_base'/' no_int_coord_base.txt
sed -i ''$lines's/d'$lines'/d'$new_labels_number_base'/' no_int_coord_base.txt
else
# Calculate new r# a# and d#
new_lines=$(($lines+1))
new_labels_number_base=$(($atoms_sugar_ion+$lines))
sed -i ''$lines's/r'$new_lines'/r'$new_labels_number_base'/' no_int_coord_base.txt
sed -i ''$lines's/a'$new_lines'/a'$new_labels_number_base'/' no_int_coord_base.txt
sed -i ''$lines's/d'$new_lines'/d'$new_labels_number_base'/' no_int_coord_base.txt
fi

# Calculate new atoms labels for the 
awk 'NR=='$lines' { print $1, $2=$2+'$atoms_sugar_ion', $3, $4=$4+'$atoms_sugar_ion', $5, $6=$6+'$atoms_sugar_ion', $7 }' no_int_coord_base.txt >> no_int_coord_base_labels_fixed.txt 

done

## Add the three first lines to the z matriz of the base and fix them accordingly
# Note: to do this is important copy first column of this file and three first rows in independent files
awk '{ print $1 }' no_int_coord_base.txt > atoms_symbols_no_int_coord_base.txt

# Print first symbol
sed -n 1,1p atoms_symbols_no_int_coord_base.txt > first_symbol.txt
first_symbol=`cat first_symbol.txt`
sed -n 2,2p atoms_symbols_no_int_coord_base.txt > second_symbol.txt
second_symbol=`cat second_symbol.txt`
sed -n 3,3p atoms_symbols_no_int_coord_base.txt > third_symbol.txt
third_symbol=`cat third_symbol.txt`
sed -n 4,4p atoms_symbols_no_int_coord_base.txt > fourth_symbol.txt
fourth_symbol=`cat fourth_symbol.txt`
sed -n 6,6p atoms_symbols_no_int_coord_base.txt > six_symbol.txt
six_symbol=`cat six_symbol.txt`
sed -n 7,7p atoms_symbols_no_int_coord_base.txt > seven_symbol.txt
seven_symbol=`cat seven_symbol.txt`
sed -n 9,9p atoms_symbols_no_int_coord_base.txt > nine_symbol.txt
nine_symbol=`cat nine_symbol.txt`


# Replace first line with internal coordinates as it should be in file
sed -i '1 s/.*/\'$first_symbol'  '$anomer'  r'$first_atom_base'  '$hemiacetal_O'  a'$first_atom_base'  '$atom_before_hemiacetal'  d'$first_atom_base'/' no_int_coord_base.txt

# Replace second line with another specific coordinates
sed -i '2 s/.*/\'$second_symbol'  '$first_atom_base'  r'$second_atom_base'  '$anomer'  a'$second_atom_base'  '$hemiacetal_O'  d'$second_atom_base'/' no_int_coord_base.txt

# Replace third line with another specific coordinates
sed -i '3 s/.*/\'$third_symbol'  '$second_atom_base'  r'$third_atom_base'  '$first_atom_base'  a'$third_atom_base'  '$anomer'  d'$third_atom_base'/' no_int_coord_base.txt

# Delete lines from line 3 in no_int_coord_base.txt
sed -i '4,$d' no_int_coord_base.txt

# Paste all together 
cat no_int_coord_base.txt no_int_coord_base_labels_fixed.txt > no_int_coord_base_labels_fixed_v2.txt

######## Replace first hydrogen line with correct atoms labels if is adenine or guanine #######

# First output number of line with first H
awk '/H/{ print NR; exit }' no_int_coord_base_labels_fixed_v2.txt > line_of_first_H.txt
line_H=`cat line_of_first_H.txt`
echo $line_H

for lines in $(seq 1 1 $ions_atoms_base)

do

if [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ] || [ "$ru" == "cytosine" ]; then

if [[ $lines -eq $line_H ]];
then
# Change the actual line
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else 
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi

elif [ "$ru" == "thymine" ] || [ "$ru" == "uracil" ]; then
if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ]; then
if [ "$lines" == 4 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi
fi

elif [ "$ru" == "CA" ]; then

if [ "$lines" == 4 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 6 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 7 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 9 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi

elif [ "$ru" == "melamine" ]; then

if [ "$lines" == 7 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 13 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi

elif [ "$ru" == "BA_cbond" ]; then
if [ "$lines" == 4 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 5 ] || [ "$lines" == 6 ] || [ "$lines" == 7 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi

elif [ "$ru" == "BA_nbond" ]; then
if [ "$lines" == 4 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 9 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 10 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi

elif [ "$ru" == "TARC_cbond" ]; then
if [ "$lines" == 5 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 8 ] || [ "$lines" == 10 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi

elif [ "$ru" == "TARC_nbond" ] && [ "$type_anomer" == "alpha" ]; then
if [ "$lines" == 4 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 6 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi

elif [ "$ru" == "TARC_nbond" ] && [ "$type_anomer" == "beta" ]; then
if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ]; then
if [ "$lines" == 4 ]; then
awk 'NR=='$lines' { print $1, $2='$hemiacetal_O', $3, $4='$atom_before_hemiacetal', $5, $6='$atom_before_before_hemiacetal', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 5 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$hemiacetal_O', $5, $6='$atom_before_hemiacetal', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 6 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 9 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi
elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ]; then
pyranose_hemiacetal_O=6
if [ "$lines" == 4 ]; then
awk 'NR=='$lines' { print $1, $2='$pyranose_hemiacetal_O', $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 5 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$pyranose_hemiacetal_O', $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 6 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 9 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$pyranose_hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi


fi

fi
done

#rm no_int_coord_base_labels_fixed_v2.txt no_int_coord_base_labels_fixed.txt no_int_coord_base.txt

# Rename this file as final
mv no_int_coord_base_labels_fixed_v3.txt no_int_coord_base_sinH.gzmat  #### THIS IS THE FINAL FILE FOR THE Z MATRIZ OF THE BASE FIXED APPROPIATELY


################ Work in internal coordinates values of base ############################

# Copy values of internal coordinates in different file
lines2=$(($ions_atoms_base+1))

# Copy lines with values of internal coordinates to another file
sed -n '/Variables:/,$p' $ru.gzmat > internal_coordinates_values_base.txt

# Delete first line of this file
sed -i '1d' internal_coordinates_values_base.txt

# Delete coordinates values for hydrogen that was deleted from base
#hydrogen_deleted=hydrogen # This have to be modify accordingly ************************(THIS LINE CAN CHANGE DEPENDING ON THE RU)***********************

sed -i '/r'$hydrogen'/d' internal_coordinates_values_base.txt
sed -i '/a'$hydrogen'/d' internal_coordinates_values_base.txt
sed -i '/d'$hydrogen'/d' internal_coordinates_values_base.txt

# Change labels of each atoms in internal coordinates file of RU before attaching to TC

for labels in $(seq 1 1 $ions_atoms_base)

do

echo "r$labels=" >> labels_base_fix.txt
echo "a$labels=" >> labels_base_fix.txt
echo "d$labels=" >> labels_base_fix.txt

done

# Eliminate everything before the = sign in old int coord file from the base
cut -f2 -d"=" internal_coordinates_values_base.txt > internal_coordinates_values_base_bef_nucleo.txt

# Paste new labels and values for internal coordinates
paste labels_base_fix.txt internal_coordinates_values_base_bef_nucleo.txt > new_internal_coordinates_base_bef_nucleo.txt

# Reduce spaces between two columns
sed -i $'s/\t//g' new_internal_coordinates_base_bef_nucleo.txt

tr -s " " < new_internal_coordinates_base_bef_nucleo.txt > internal_coordinates_base_sinH_before.gzmat
# Delete previous file
rm internal_coordinates_values_base.txt

cp internal_coordinates_base_sinH_before.gzmat internal_coordinates_values_base.txt

# Add new lines with new values of distances, angles and dihedrals to file with internal coordinates (this is where I add the numbers of internal coordinates to keep in all molecules)*****IMPORTANT!

# This one allways remains the same

if [ "$type_anomer" == "alpha" ];

then

if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ];

then

# Insert a1 and d1 for purines
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 98.26" >> first_int_coord_base.txt
echo "d1= 99.56" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for purines
sed -i '2 i a2= 132.67' internal_coordinates_values_base.txt
sed -i '3 i d2= 216.01' internal_coordinates_values_base.txt
sed -i '6 i d3= 185.94' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$sugar_ring" == "D1C4_alpha" ] && [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ]; then

# Insert a1 and d1 for purines
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 116.43" >> first_int_coord_base.txt
echo "d1= 179.29" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for purines
sed -i '2 i a2= 129.09' internal_coordinates_values_base.txt
sed -i '3 i d2= 169.90' internal_coordinates_values_base.txt
sed -i '6 i d3= 176.15' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$sugar_ring" == "D4C1_alpha" ] && [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ]; then

# Insert a1 and d1 for purines
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 112.99" >> first_int_coord_base.txt
echo "d1= 293.05" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for purines
sed -i '2 i a2= 128.06' internal_coordinates_values_base.txt
sed -i '3 i d2= 165.71' internal_coordinates_values_base.txt
sed -i '6 i d3= 179.29' internal_coordinates_values_base.txt

elif [ "$ru" == "cytosine" ] || [ "$ru" == "thymine" ] || [ "$ru" == "uracil" ]; then

if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ]; then
# Insert a1 and d1 for pyrimidines
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 111.67" >> first_int_coord_base.txt
echo "d1=  97.87" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 109.63' internal_coordinates_values_base.txt
sed -i '3 i d2=  45.02' internal_coordinates_values_base.txt
sed -i '6 i d3= 354.43' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$sugar_ring" == "D1C4_alpha" ]; then

# Insert a1 and d1 for pyrimidines
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 114.98" >> first_int_coord_base.txt
echo "d1= 178.69" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 115.03' internal_coordinates_values_base.txt
sed -i '3 i d2= 168.42' internal_coordinates_values_base.txt
sed -i '6 i d3= 357.60' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$sugar_ring" == "D4C1_alpha" ]; then

# Insert a1 and d1 for pyrimidines
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 111.29" >> first_int_coord_base.txt
echo "d1= 292.09" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 115.53' internal_coordinates_values_base.txt
sed -i '3 i d2= 165.94' internal_coordinates_values_base.txt
sed -i '6 i d3= 359.84' internal_coordinates_values_base.txt

fi

elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "CA" ];

then

# Insert a, d for non canonical base
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 111.67" >> first_int_coord_base.txt # a1=113.40
echo "d1=  97.87" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 109.63' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2=  45.02' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 354.43' internal_coordinates_values_base.txt # d3=0.40

sed -i '8s/.*/a4= 113.17/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 180.40/' internal_coordinates_values_base.txt

sed -i '14s/.*/a6= 116.65/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6= 223.45/' internal_coordinates_values_base.txt

sed -i '17s/.*/a7= 139.29/' internal_coordinates_values_base.txt
sed -i '18s/.*/d7= 179.60/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 123.42/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9= 359.59/' internal_coordinates_values_base.txt


elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$ru" == "CA" ]; then

if [ "$sugar_ring" == "D1C4_alpha" ]; then

# Insert a, d for non canonical base
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 115.40" >> first_int_coord_base.txt # a1=113.40
echo "d1= 178.29" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 116.53' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 168.12' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 357.87' internal_coordinates_values_base.txt # d3=0.40

sed -i '8s/.*/a4= 113.17/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 177.87/' internal_coordinates_values_base.txt

sed -i '14s/.*/a6= 116.59/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6= 346.22/' internal_coordinates_values_base.txt

sed -i '17s/.*/a7= 139.28/' internal_coordinates_values_base.txt
sed -i '18s/.*/d7= 182.12/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 123.41/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9=   2.13/' internal_coordinates_values_base.txt

elif [ "$sugar_ring" == "D4C1_alpha" ]; then

# Insert a, d for non canonical base
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 111.25" >> first_int_coord_base.txt # a1=113.40
echo "d1= 292.12" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 115.81' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 164.73' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3=   1.16' internal_coordinates_values_base.txt # d3=0.40

sed -i '8s/.*/a4= 113.18/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 181.16/' internal_coordinates_values_base.txt

sed -i '14s/.*/a6= 117.34/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6= 345.77/' internal_coordinates_values_base.txt

sed -i '17s/.*/a7= 139.29/' internal_coordinates_values_base.txt
sed -i '18s/.*/d7= 178.82/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 123.41/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9= 358.83/' internal_coordinates_values_base.txt

fi

elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "melamine" ];

then

# Insert a, d for non canonical base
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 112.66" >> first_int_coord_base.txt # a1=113.40
echo "d1=  89.85" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 117.45' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2=  18.74' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 191.82' internal_coordinates_values_base.txt # d3=0.40

sed -i '18s/.*/d7=  12.91/' internal_coordinates_values_base.txt

sed -i '36s/.*/d13= 223.16/' internal_coordinates_values_base.txt


elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$ru" == "melamine" ]; then

if [ "$sugar_ring" == "D1C4_alpha" ]; then

# Insert a, d for non canonical base
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 115.20" >> first_int_coord_base.txt # a1=113.40
echo "d1= 178.19" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 116.43' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 143.07' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 190.30' internal_coordinates_values_base.txt # d3=0.40

sed -i '18s/.*/d7=  11.39/' internal_coordinates_values_base.txt

sed -i '36s/.*/d13= 346.32/' internal_coordinates_values_base.txt

elif [ "$sugar_ring" == "D4C1_alpha" ]; then

# Insert a, d for non canonical base
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 111.00" >> first_int_coord_base.txt # a1=113.40
echo "d1= 291.87" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 117.99' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 141.69' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 191.79' internal_coordinates_values_base.txt # d3=0.40

sed -i '18s/.*/d7= 359.22/' internal_coordinates_values_base.txt

sed -i '36s/.*/d13= 345.97/' internal_coordinates_values_base.txt

fi

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$ru" == "BA_cbond" ]; then

if [ "$sugar_ring" == "D1C4_alpha" ]; then

# Insert a, d for BA_cbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 115.44" >> first_int_coord_base.txt # a1=113.40
echo "d1= 178.48" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_cbond 
sed -i '2 i a2= 117.80' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 346.10' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 358.96' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for BA_cbond
sed -i '8s/.*/a4= 122.60/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 165.02/' internal_coordinates_values_base.txt

sed -i '12s/.*/d5= 181.10/' internal_coordinates_values_base.txt

sed -i '15s/.*/d6= 181.10/' internal_coordinates_values_base.txt

sed -i '17s/.*/a7= 126.57/' internal_coordinates_values_base.txt
sed -i '18s/.*/d7=   1.10/' internal_coordinates_values_base.txt

elif [ "$sugar_ring" == "D4C1_alpha" ]; then

# Insert a, d for BA_cbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 110.66" >> first_int_coord_base.txt # a1=113.40
echo "d1= 291.71" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_cbond 
sed -i '2 i a2= 118.38' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 346.13' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 359.83' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for BA_cbond
sed -i '8s/.*/a4= 122.60/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 165.02/' internal_coordinates_values_base.txt

sed -i '12s/.*/d5= 180.18/' internal_coordinates_values_base.txt

sed -i '15s/.*/d6= 180.18/' internal_coordinates_values_base.txt

sed -i '17s/.*/a7= 126.57/' internal_coordinates_values_base.txt
sed -i '18s/.*/d7=   0.18/' internal_coordinates_values_base.txt

fi

elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "BA_cbond" ]; then

# Insert a, d for BA_cbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 114.06" >> first_int_coord_base.txt # a1=113.40
echo "d1=  89.71" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_cbond 
sed -i '2 i a2= 119.03' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 223.63' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3=   0.40' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for BA_cbond
sed -i '8s/.*/a4= 121.37/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4=  44.04/' internal_coordinates_values_base.txt

sed -i '12s/.*/d5= 179.59/' internal_coordinates_values_base.txt

sed -i '15s/.*/d6= 179.58/' internal_coordinates_values_base.txt

sed -i '17s/.*/a7= 126.56/' internal_coordinates_values_base.txt
sed -i '18s/.*/d7= 359.59/' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "BA_nbond" ]; then

# Insert a, d for BA_nbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 112.98" >> first_int_coord_base.txt 
echo "d1=  89.36" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for BA_nbond 
sed -i '2 i a2= 116.02' internal_coordinates_values_base.txt 
sed -i '3 i d2= 223.53' internal_coordinates_values_base.txt 
sed -i '6 i d3= 359.47' internal_coordinates_values_base.txt 

# Insert other dihedral angles for BA_nbond
sed -i '8s/.*/a4= 84.97/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 179.46/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 116.22/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9=  43.05/' internal_coordinates_values_base.txt

sed -i '26s/.*/a10= 124.78/' internal_coordinates_values_base.txt
sed -i '27s/.*/d10=   0.54/' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$ru" == "BA_nbond" ]; then

if [ "$sugar_ring" == "D1C4_alpha" ]; then

# Insert a, d for BA_nbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 115.27" >> first_int_coord_base.txt 
echo "d1= 178.65" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for BA_nbond 
sed -i '2 i a2= 116.10' internal_coordinates_values_base.txt 
sed -i '3 i d2= 346.05' internal_coordinates_values_base.txt 
sed -i '6 i d3=   0.27' internal_coordinates_values_base.txt 

# Insert other dihedral angles for BA_nbond
sed -i '8s/.*/a4=  84.97/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 180.26/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 116.13/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9= 166.28/' internal_coordinates_values_base.txt

sed -i '26s/.*/a10= 124.78/' internal_coordinates_values_base.txt
sed -i '27s/.*/d10= 359.73/' internal_coordinates_values_base.txt

elif [ "$sugar_ring" == "D4C1_alpha" ]; then

# Insert a, d for BA_nbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 111.03" >> first_int_coord_base.txt 
echo "d1= 291.96" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for BA_nbond 
sed -i '2 i a2= 116.34' internal_coordinates_values_base.txt 
sed -i '3 i d2= 345.91' internal_coordinates_values_base.txt 
sed -i '6 i d3= 359.53' internal_coordinates_values_base.txt 

# Insert other dihedral angles for BA_nbond
sed -i '8s/.*/a4=  84.97/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 179.54/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 115.90/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9= 165.51/' internal_coordinates_values_base.txt

sed -i '26s/.*/a10= 124.78/' internal_coordinates_values_base.txt
sed -i '27s/.*/d10=   0.46/' internal_coordinates_values_base.txt

fi

elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "TARC_cbond" ]; then
# Insert a, d for TARC_cbond
echo "r1= 1.5199" > first_int_coord_base.txt
echo "a1= 113.46" >> first_int_coord_base.txt 
echo "d1=  89.02" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for TARC_cbond 
sed -i '2 i a2= 149.52' internal_coordinates_values_base.txt  
sed -i '3 i d2= 181.28' internal_coordinates_values_base.txt  
sed -i '6 i d3= 179.63' internal_coordinates_values_base.txt 

# Insert other dihedral angles for TARC_cbond
sed -i '11s/.*/a5= 122.48/' internal_coordinates_values_base.txt
sed -i '12s/.*/d5=   1.62/' internal_coordinates_values_base.txt

sed -i '21s/.*/d8=   0.24/' internal_coordinates_values_base.txt

sed -i '27s/.*/d10=  2.30/' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$ru" == "TARC_cbond" ]; then

if [ "$sugar_ring" == "D1C4_alpha" ]; then

# Insert a, d for TARC_cbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 113.74" >> first_int_coord_base.txt 
echo "d1= 177.87" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for TARC_cbond 
sed -i '2 i a2= 150.20' internal_coordinates_values_base.txt  
sed -i '3 i d2= 345.90' internal_coordinates_values_base.txt  
sed -i '6 i d3= 180.99' internal_coordinates_values_base.txt 

# Insert other dihedral angles for TARC_cbond
sed -i '11s/.*/a5= 121.78/' internal_coordinates_values_base.txt
sed -i '12s/.*/d5= 167.84/' internal_coordinates_values_base.txt

sed -i '21s/.*/d8=   1.60/' internal_coordinates_values_base.txt

sed -i '27s/.*/d10=  1.50/' internal_coordinates_values_base.txt

elif [ "$sugar_ring" == "D4C1_alpha" ]; then

# Insert a, d for TARC_cbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 111.61" >> first_int_coord_base.txt 
echo "d1= 292.35" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for TARC_cbond 
sed -i '2 i a2= 149.38' internal_coordinates_values_base.txt  
sed -i '3 i d2= 303.77' internal_coordinates_values_base.txt  
sed -i '6 i d3= 178.29' internal_coordinates_values_base.txt 

# Insert other dihedral angles for TARC_cbond
sed -i '11s/.*/a5= 122.61/' internal_coordinates_values_base.txt
sed -i '12s/.*/d5= 122.51/' internal_coordinates_values_base.txt

sed -i '21s/.*/d8= 358.89/' internal_coordinates_values_base.txt

sed -i '27s/.*/d10=  3.11/' internal_coordinates_values_base.txt

fi

elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "TARC_nbond" ]; then
# Insert a, d for TARC_nbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 112.76" >> first_int_coord_base.txt 
echo "d1=  89.85" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for TARC_nbond   
sed -i '2 i a2= 118.80' internal_coordinates_values_base.txt  
sed -i '3 i d2=  77.95' internal_coordinates_values_base.txt  
sed -i '6 i d3= 160.79' internal_coordinates_values_base.txt 

sed -i '9s/.*/d4= 339.95/' internal_coordinates_values_base.txt 

sed -i '14s/.*/a6= 116.47/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6= 223.15/' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$ru" == "TARC_nbond" ]; then
if [ "$sugar_ring" == "D1C4_alpha" ]; then

# Insert a, d for TARC_nbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 113.49" >> first_int_coord_base.txt 
echo "d1= 177.50" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for TARC_nbond   
sed -i '2 i a2= 116.75' internal_coordinates_values_base.txt  
sed -i '3 i d2= 205.73' internal_coordinates_values_base.txt  
sed -i '6 i d3= 156.72' internal_coordinates_values_base.txt 

sed -i '9s/.*/d4= 335.88/' internal_coordinates_values_base.txt 

sed -i '14s/.*/a6= 115.95/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6= 347.03/' internal_coordinates_values_base.txt

elif [ "$sugar_ring" == "D4C1_alpha" ]; then

# Insert a, d for TARC_nbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 110.57" >> first_int_coord_base.txt 
echo "d1= 290.85" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for TARC_nbond   
sed -i '2 i a2= 118.78' internal_coordinates_values_base.txt  
sed -i '3 i d2= 204.45' internal_coordinates_values_base.txt  
sed -i '6 i d3= 157.01' internal_coordinates_values_base.txt 

sed -i '9s/.*/d4= 336.17/' internal_coordinates_values_base.txt 

sed -i '14s/.*/a6= 114.51/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6= 346.57/' internal_coordinates_values_base.txt

fi

fi






else # This is where the values of the r, a, d in internal coord of bases for beta sugars are changed






if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ];

then

# Insert a1 and d1 for purines
echo "r1=   1.52" > first_int_coord_base.txt
echo "a1= 117.78" >> first_int_coord_base.txt
echo "d1= 223.12" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for purines
sed -i '2 i a2= 131.20' internal_coordinates_values_base.txt
sed -i '3 i d2= 334.71' internal_coordinates_values_base.txt
sed -i '6 i d3= 177.00' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$sugar_ring" == "D1C4_beta" ] && [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ]; then

# Insert a1 and d1 for purines
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 110.09" >> first_int_coord_base.txt
echo "d1= 68.20" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for purines
sed -i '2 i a2= 127.70' internal_coordinates_values_base.txt
sed -i '3 i d2= 227.82' internal_coordinates_values_base.txt
sed -i '6 i d3= 180.34' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$sugar_ring" == "D4C1_beta" ] && [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ]; then

# Insert a1 and d1 for purines
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 116.85" >> first_int_coord_base.txt
echo "d1= 178.85" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for purines
sed -i '2 i a2= 128.36' internal_coordinates_values_base.txt
sed -i '3 i d2= 229.04' internal_coordinates_values_base.txt
sed -i '6 i d3= 174.76' internal_coordinates_values_base.txt

elif [ "$ru" == "cytosine" ] || [ "$ru" == "thymine" ] || [ "$ru" == "uracil" ]; then

if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ]; then

# Insert a1 and d1 for pyrimidines
echo "r1=   1.52" > first_int_coord_base.txt
echo "a1= 108.07" >> first_int_coord_base.txt
echo "d1= 201.62" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 116.84' internal_coordinates_values_base.txt
sed -i '3 i d2= 170.38' internal_coordinates_values_base.txt
sed -i '6 i d3= 359.48' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$sugar_ring" == "D1C4_beta" ]; then

# Insert a1 and d1 for pyrimidines
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 109.60" >> first_int_coord_base.txt
echo "d1=  65.99" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 115.90' internal_coordinates_values_base.txt
sed -i '3 i d2= 228.78' internal_coordinates_values_base.txt
sed -i '6 i d3= 359.19' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$sugar_ring" == "D4C1_beta" ]; then

# Insert a1 and d1 for pyrimidines
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 114.60" >> first_int_coord_base.txt
echo "d1= 180.21" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 115.64' internal_coordinates_values_base.txt
sed -i '3 i d2= 227.79' internal_coordinates_values_base.txt
sed -i '6 i d3= 356.80' internal_coordinates_values_base.txt

fi

elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "CA" ];

then

# Insert a, d for non canonical base
echo "r1=   1.52" > first_int_coord_base.txt
echo "a1= 111.67" >> first_int_coord_base.txt # a1=113.40
echo "d1= 202.71" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 114.33' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 348.91' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 3.40' internal_coordinates_values_base.txt # d3=0.40

sed -i '8s/.*/a4= 113.17/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 183.39/' internal_coordinates_values_base.txt

sed -i '14s/.*/a6= 118.75/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6= 172.06/' internal_coordinates_values_base.txt

sed -i '17s/.*/a7= 139.29/' internal_coordinates_values_base.txt
sed -i '18s/.*/d7= 176.47/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 123.42/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9= 356.45/' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ]  && [ "$ru" == "CA" ]; then

if [ "$sugar_ring" == "D1C4_beta" ]; then

# Insert a, d for non canonical base
echo "r1=   1.52" > first_int_coord_base.txt
echo "a1= 110.69" >> first_int_coord_base.txt # a1=113.40
echo "d1=  67.46" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 118.10' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 228.34' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 359.81' internal_coordinates_values_base.txt # d3=0.40

sed -i '8s/.*/a4= 113.17/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 179.81/' internal_coordinates_values_base.txt

sed -i '14s/.*/a6= 115.06/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6=  48.18/' internal_coordinates_values_base.txt

sed -i '17s/.*/a7= 139.29/' internal_coordinates_values_base.txt
sed -i '18s/.*/d7= 180.19/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 123.41/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9=   0.18/' internal_coordinates_values_base.txt

elif [ "$sugar_ring" == "D4C1_beta" ]; then

# Insert a, d for non canonical base
echo "r1=   1.52" > first_int_coord_base.txt
echo "a1= 115.33" >> first_int_coord_base.txt # a1=113.40
echo "d1= 179.97" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 115.45' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 228.11' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3=   4.03' internal_coordinates_values_base.txt # d3=0.40

sed -i '8s/.*/a4= 113.18/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 184.03/' internal_coordinates_values_base.txt

sed -i '14s/.*/a6= 117.59/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6=  51.74/' internal_coordinates_values_base.txt

sed -i '17s/.*/a7= 139.29/' internal_coordinates_values_base.txt
sed -i '18s/.*/d7= 175.90/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 123.42/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9= 355.90/' internal_coordinates_values_base.txt

fi

elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "melamine" ];

then

# Insert a, d for non canonical base
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 110.96" >> first_int_coord_base.txt # a1=113.40
echo "d1= 203.29" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 121.79' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 144.42' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 192.65' internal_coordinates_values_base.txt # d3=0.40

sed -i '18s/.*/d7=  13.75/' internal_coordinates_values_base.txt

sed -i '36s/.*/d13= 348.66/' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$ru" == "melamine" ]; then

if [ "$sugar_ring" == "D1C4_beta" ]; then

# Insert a, d for non canonical base
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 110.96" >> first_int_coord_base.txt # a1=113.40
echo "d1=  67.67" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 118.25' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2=  24.94' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 190.89' internal_coordinates_values_base.txt # d3=0.40

sed -i '18s/.*/d7=  11.98/' internal_coordinates_values_base.txt

sed -i '36s/.*/d13= 228.31/' internal_coordinates_values_base.txt

elif [ "$sugar_ring" == "D4C1_beta" ]; then

# Insert a, d for non canonical base
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 115.08" >> first_int_coord_base.txt # a1=113.40
echo "d1= 179.79" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 117.65' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2=  20.56' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 195.34' internal_coordinates_values_base.txt # d3=0.40

sed -i '18s/.*/d7=  16.44/' internal_coordinates_values_base.txt

sed -i '36s/.*/d13= 228.12/' internal_coordinates_values_base.txt

fi

elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "BA_cbond" ]; then

# Insert a, d for BA_cbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 109.50" >> first_int_coord_base.txt # a1=113.40
echo "d1= 202.57" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_cbond 
sed -i '2 i a2= 117.28' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 167.34' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3=   1.77' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for BA_cbond
sed -i '8s/.*/a4= 123.10/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 349.18/' internal_coordinates_values_base.txt

sed -i '12s/.*/d5= 178.12/' internal_coordinates_values_base.txt

sed -i '15s/.*/d6= 178.12/' internal_coordinates_values_base.txt

sed -i '18s/.*/d7= 358.12/' internal_coordinates_values_base.txt

sed -i '24s/.*/d9= 359.97/' internal_coordinates_values_base.txt

sed -i '30s/.*/d11= 179.97/' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$ru" == "BA_cbond" ]; then

if [ "$sugar_ring" == "D1C4_beta" ]; then

# Insert a, d for BA_cbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 110.81" >> first_int_coord_base.txt # a1=113.40
echo "d1=  67.41" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_cbond 
sed -i '2 i a2= 117.16' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 228.40' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3=   1.54' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for BA_cbond
sed -i '8s/.*/a4= 123.23/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4=  50.00/' internal_coordinates_values_base.txt

sed -i '12s/.*/d5= 178.36/' internal_coordinates_values_base.txt

sed -i '15s/.*/d6= 178.36/' internal_coordinates_values_base.txt

sed -i '17s/.*/a7= 126.56/' internal_coordinates_values_base.txt
sed -i '18s/.*/d7= 358.36/' internal_coordinates_values_base.txt

elif [ "$sugar_ring" == "D4C1_beta" ]; then

# Insert a, d for BA_cbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 113.97" >> first_int_coord_base.txt # a1=113.40
echo "d1= 180.84" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_cbond 
sed -i '2 i a2= 117.64' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 227.33' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3=   0.91' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for BA_cbond
sed -i '8s/.*/a4= 122.76/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4=  48.27/' internal_coordinates_values_base.txt

sed -i '12s/.*/d5= 179.04/' internal_coordinates_values_base.txt

sed -i '15s/.*/d6= 179.04/' internal_coordinates_values_base.txt

sed -i '17s/.*/a7= 126.57/' internal_coordinates_values_base.txt
sed -i '18s/.*/d7= 359.04/' internal_coordinates_values_base.txt

fi


elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "BA_nbond" ]; then

# Insert a, d for BA_nbond
echo "r1= 1.5201" > first_int_coord_base.txt
echo "a1= 109.12" >> first_int_coord_base.txt # a1=113.40
echo "d1= 202.19" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_nbond 
sed -i '2 i a2= 117.00' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 349.42' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 357.46' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for BA_nbond
sed -i '8s/.*/a4=  84.97/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 177.46/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 115.20/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9= 167.20/' internal_coordinates_values_base.txt

sed -i '27s/.*/d10=  2.51/' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$ru" == "BA_nbond" ]; then

if [ "$sugar_ring" == "D1C4_beta" ]; then

# Insert a, d for BA_nbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 110.80" >> first_int_coord_base.txt 
echo "d1=  68.01" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for BA_nbond 
sed -i '2 i a2= 115.48' internal_coordinates_values_base.txt 
sed -i '3 i d2= 228.09' internal_coordinates_values_base.txt 
sed -i '6 i d3=   0.06' internal_coordinates_values_base.txt 

# Insert other dihedral angles for BA_nbond
sed -i '8s/.*/a4=  84.97/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 180.03/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 116.76/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9=  48.14/' internal_coordinates_values_base.txt

sed -i '26s/.*/a10= 124.78/' internal_coordinates_values_base.txt
sed -i '27s/.*/d10= 359.94/' internal_coordinates_values_base.txt

elif [ "$sugar_ring" == "D4C1_beta" ]; then

# Insert a, d for BA_nbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 114.57" >> first_int_coord_base.txt 
echo "d1= 180.69" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for BA_nbond 
sed -i '2 i a2= 114.75' internal_coordinates_values_base.txt 
sed -i '3 i d2= 227.54' internal_coordinates_values_base.txt 
sed -i '6 i d3=   0.92' internal_coordinates_values_base.txt 

# Insert other dihedral angles for BA_nbond
sed -i '8s/.*/a4=  84.97/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 180.93/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 117.48/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9=  48.37/' internal_coordinates_values_base.txt

sed -i '26s/.*/a10= 124.78/' internal_coordinates_values_base.txt
sed -i '27s/.*/d10= 359.04/' internal_coordinates_values_base.txt

fi

elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "TARC_cbond" ]; then
# Insert a, d for TARC_cbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 110.52" >> first_int_coord_base.txt 
echo "d1= 203.37" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for TARC_cbond 
sed -i '2 i a2= 149.69' internal_coordinates_values_base.txt  
sed -i '3 i d2= 168.33' internal_coordinates_values_base.txt  
sed -i '6 i d3= 179.54' internal_coordinates_values_base.txt 

# Insert other dihedral angles for TARC_cbond
sed -i '11s/.*/a5= 122.31/' internal_coordinates_values_base.txt
sed -i '12s/.*/d5= 348.58/' internal_coordinates_values_base.txt

# Insert other dihedral angles for TARC_cbond
sed -i '21s/.*/d8=   0.16/' internal_coordinates_values_base.txt

sed -i '27s/.*/d10=  2.34/' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$ru" == "TARC_cbond" ]; then

if [ "$sugar_ring" == "D1C4_beta" ]; then

# Insert a, d for TARC_cbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 108.37" >> first_int_coord_base.txt 
echo "d1=  64.74" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for TARC_cbond 
sed -i '2 i a2= 150.18' internal_coordinates_values_base.txt  
sed -i '3 i d2=  50.36' internal_coordinates_values_base.txt  
sed -i '6 i d3= 178.22' internal_coordinates_values_base.txt 

# Insert other dihedral angles for TARC_cbond
sed -i '11s/.*/a5= 121.81/' internal_coordinates_values_base.txt
sed -i '12s/.*/d5= 229.04/' internal_coordinates_values_base.txt

sed -i '21s/.*/d8= 358.82/' internal_coordinates_values_base.txt

sed -i '27s/.*/d10=  3.13/' internal_coordinates_values_base.txt

elif [ "$sugar_ring" == "D4C1_beta" ]; then

# Insert a, d for TARC_cbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 114.33" >> first_int_coord_base.txt 
echo "d1= 179.97" >> first_int_coord_base.txt 

# Insert a2 d2 and d3 for TARC_cbond 
sed -i '2 i a2= 151.02' internal_coordinates_values_base.txt  
sed -i '3 i d2=  52.24' internal_coordinates_values_base.txt  
sed -i '6 i d3= 175.60' internal_coordinates_values_base.txt 

# Insert other dihedral angles for TARC_cbond
sed -i '11s/.*/a5= 120.90/' internal_coordinates_values_base.txt
sed -i '12s/.*/d5= 227.88/' internal_coordinates_values_base.txt

sed -i '21s/.*/d8= 356.21/' internal_coordinates_values_base.txt

sed -i '27s/.*/d10=  4.57/' internal_coordinates_values_base.txt

fi

elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "TARC_nbond" ] && [ "$sugar_ring" == "2endo_beta" ]; then

# Insert a, d for BA_nbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 112.00" >> first_int_coord_base.txt # a1=113.40
echo "d1= 203.67" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_nbond 
sed -i '2 i a2= 118.71' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 348.13' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 185.56' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for TARC_cbond
sed -i '7s/.*/r4= 1.7905/' internal_coordinates_values_base.txt
sed -i '8s/.*/a4= 142.29/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 280.77/' internal_coordinates_values_base.txt

# Insert other dihedral angles for TARC_cbond
sed -i '11s/.*/a5= 114.09/' internal_coordinates_values_base.txt
sed -i '12s/.*/d5= 147.57/' internal_coordinates_values_base.txt

# Insert other dihedral angles for TARC_cbond
sed -i '14s/.*/a6= 117.55/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6= 135.06/' internal_coordinates_values_base.txt

sed -i '24s/.*/d9= 216.53/' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] && [ "$ru" == "TARC_nbond" ] && [ "$sugar_ring" == "3endo_beta" ]; then

# Insert a, d for BA_nbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 109.76" >> first_int_coord_base.txt # a1=113.40
echo "d1= 222.29" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_nbond 
sed -i '2 i a2= 118.19' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 350.54' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 156.42' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for TARC_cbond
sed -i '7s/.*/r4= 1.8440/' internal_coordinates_values_base.txt
sed -i '8s/.*/a4= 148.69/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 335.55/' internal_coordinates_values_base.txt

# Insert other dihedral angles for TARC_cbond
sed -i '11s/.*/a5= 108.28/' internal_coordinates_values_base.txt
sed -i '12s/.*/d5= 116.43/' internal_coordinates_values_base.txt

# Insert other dihedral angles for TARC_cbond
sed -i '14s/.*/a6= 114.66/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6= 132.03/' internal_coordinates_values_base.txt

sed -i '24s/.*/d9= 221.27/' internal_coordinates_values_base.txt

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] && [ "$ru" == "TARC_nbond" ]; then

if [ "$sugar_ring" == "D1C4_beta" ]; then

# Insert a, d for BA_nbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 111.18" >> first_int_coord_base.txt # a1=113.40
echo "d1=  68.16" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_nbond 
sed -i '2 i a2= 118.60' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 228.11' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 160.71' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for TARC_cbond
sed -i '7s/.*/r4= 1.7970/' internal_coordinates_values_base.txt
sed -i '8s/.*/a4= 105.35/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 161.53/' internal_coordinates_values_base.txt

# Insert other dihedral angles for TARC_cbond
sed -i '11s/.*/a5= 115.26/' internal_coordinates_values_base.txt
sed -i '12s/.*/d5= 303.69/' internal_coordinates_values_base.txt

# Insert other dihedral angles for TARC_cbond
sed -i '14s/.*/a6= 116.59/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6=  13.20/' internal_coordinates_values_base.txt

sed -i '24s/.*/d9= 215.56/' internal_coordinates_values_base.txt

elif [ "$sugar_ring" == "D4C1_beta" ]; then

# Insert a, d for BA_nbond
echo "r1= 1.5200" > first_int_coord_base.txt
echo "a1= 114.10" >> first_int_coord_base.txt # a1=113.40
echo "d1= 180.74" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_nbond 
sed -i '2 i a2= 118.28' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 227.41' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 159.88' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for TARC_cbond
sed -i '7s/.*/r4= 1.7396/' internal_coordinates_values_base.txt
sed -i '8s/.*/a4= 108.23/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 169.87/' internal_coordinates_values_base.txt

# Insert other dihedral angles for TARC_cbond
sed -i '11s/.*/a5= 109.92/' internal_coordinates_values_base.txt
sed -i '12s/.*/d5= 293.70/' internal_coordinates_values_base.txt

# Insert other dihedral angles for TARC_cbond
sed -i '14s/.*/a6= 116.43/' internal_coordinates_values_base.txt
sed -i '15s/.*/d6=  11.72/' internal_coordinates_values_base.txt

sed -i '24s/.*/d9= 218.42/' internal_coordinates_values_base.txt

fi

fi

fi

# Put both files together
cat first_int_coord_base.txt internal_coordinates_values_base.txt > internal_coordinates_base.txt

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

# Fix value of dihedral angle between first hydrogen atom of base:

# First print dihedral angle appropiately
awk '(NR=='$line_H'){print $7}' no_int_coord_base_sinH.gzmat > dihedral_level_firstH.txt
value=`cat dihedral_level_firstH.txt`

# Output number of line in internal coordinate file with this value
grep -n "$value" internal_coordinates_base_sinH.gzmat | cut -d : -f1 > line_number_in_intcoord_base.txt

line=`cat line_number_in_intcoord_base.txt`

# Replace the line
if [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ] || [ "$ru" == "cytosine" ] && [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ]; then
awk 'NR=='$line' {$0="'$value'= 4.26"} 1' internal_coordinates_base_sinH.gzmat > internal_coordinates_base_sinH_v2.gzmat
elif [ "$ru" == "thymine" ] || [ "$ru" == "uracil" ] && [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "threose" ] || [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ]; then
awk 'NR=='$line' {$0="'$value'= 181.26"} 1' internal_coordinates_base_sinH.gzmat > internal_coordinates_base_sinH_v2.gzmat
elif [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ] || [ "$ru" == "cytosine" ] && [ "$sugar_ring" == "D1C4_alpha" ] || [ "$sugar_ring" == "D1C4_beta" ] && [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ]; then
awk 'NR=='$line' {$0="'$value'= 356.15"} 1' internal_coordinates_base_sinH.gzmat > internal_coordinates_base_sinH_v2.gzmat
elif [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ] || [ "$ru" == "cytosine" ] && [ "$sugar_ring" == "D4C1_alpha" ] || [ "$sugar_ring" == "D4C1_beta" ] && [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ]; then
awk 'NR=='$line' {$0="'$value'= 359.29"} 1' internal_coordinates_base_sinH.gzmat > internal_coordinates_base_sinH_v2.gzmat


fi


################### Now paste all together to create final z matriz file for the nucleosides ######################
if [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ] || [ "$ru" == "cytosine" ] || [ "$ru" == "thymine" ] || [ "$ru" == "uracil" ]; then
cat no_int_coord_sugar_sinOH.gzmat no_int_coord_base_sinH.gzmat internal_coordinates_values_sugar_sinOH.txt internal_coordinates_base_sinH_v2.gzmat > nucleoside.gzmat
elif [ "$ru" == "CA" ] || [ "$ru" == "melamine" ] || [ "$ru" == "BA_cbond" ] || [ "$ru" == "BA_nbond" ] || [ "$ru" == "TARC_cbond" ] || [ "$ru" == "TARC_nbond" ]; then
cat no_int_coord_sugar_sinOH.gzmat no_int_coord_base_sinH.gzmat internal_coordinates_values_sugar_sinOH.txt internal_coordinates_base_sinH.gzmat > nucleoside.gzmat
fi

# Create new name of nucleoside
#read -p "insert name of your output file gjf:" nucleoside
if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "2deoxy_ribopyranose" ]; then
nucleoside=2deoxy_"$ru"_"$sugar_ring"
elif [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "ribopyranose" ]; then
nucleoside=ribo_"$ru"_"$sugar_ring"
elif [ "$type_sugar" == "threose" ]; then
nucleoside=threo_"$ru"_"$sugar_ring"
fi


# Convert this file to gjf
obabel nucleoside.gzmat -O $nucleoside.gjf

# Paste instructions for b3lyp
sed -i -e 1,2d $nucleoside.gjf

cat instructions_DFT/instructions_b3lyp.txt $nucleoside.gjf > "$nucleoside"_new.gjf

# Replace name chk by name of cell
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$nucleoside'.chk/g' "$nucleoside"_new.gjf

# Add frozen coordinates
if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ]; then
cat "$nucleoside"_new.gjf frozen_coordinates/frozen_coord_fura.txt > "$nucleoside"_new2.gjf
elif [ "$type_sugar" == "threose" ]; then
cat "$nucleoside"_new.gjf frozen_coordinates/frozen_coord_threo.txt > "$nucleoside"_new2.gjf
elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ]; then
cat "$nucleoside"_new.gjf frozen_coordinates/frozen_coord_pyra.txt > "$nucleoside"_new2.gjf
fi



rm $nucleoside.gjf "$nucleoside"_new.gjf

mv "$nucleoside"_new2.gjf $nucleoside.gjf

#rm *.txt
#rm *.log


