#!/bin/bash

# Note1: no_int_coord_base.gzmat and no_int_coord_sugar.gzmat are the files with only the z matrix
# Note2: this initial test was with beta 2deoxyribofuranose and guanine

# Insert the names of the files for the TC and the RUs
read -p "insert name of log file for the Trifunctional Connector (TC) without extension:" tc
read -p "insert name of log file for the Recognition unit (RU) without extension:" ru
read -p "insert type of anomer in case that it is (e.g. beta, or alpha, or none if there is not anomer)" type_anomer
read -p "insert type of sugar in case that it is (e.g. 2deoxy_ribofuranose, or ribofuranose, or 2deoxy_ribopyranose, or pyranose, or threose, or none if the TC is not a sugar)" type_sugar


# Convert both files for the sugar and base to gzmat
obabel *.log -O *.gzmat

########################
## WORK ON SUGAR FILE ##
########################

# Declare labels of anomeric carbon
read -p "insert label of anomeric carbon as appears in gaussview:" anomer ####### ************************(THIS LINE CAN CHANGE DEPENDING ON THE TC)***********************

# Declare label of hemiacetalic oxygen for beta or 
read -p "insert label of hemiacetalic oxygen as appears in gaussview:" hemiacetal_O ####### ************************(THIS LINE CAN CHANGE DEPENDING ON THE TC)***********************

atom_before_hemiacetal=$(($hemiacetal_O-1))

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

awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt

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

elif [ "$ru" == "melamine" ] && [ "$type_anomer" == "alpha" ]; then

if [ "$lines" == 5 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 9 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 12 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 13 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 14 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 8 ]; then
awk 'NR=='$lines' { print $1, $2=$2-2, $3, $4=$4+1, $5, $6=$6-1, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi

elif [ "$ru" == "melamine" ] && [ "$type_anomer" == "beta" ]; then
if [ "$lines" == 12 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 13 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 14 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 5 ]; then
awk 'NR=='$lines' { print $1, $2=$2+2, $3, $4, $5, $6=$6-2, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 8 ]; then
awk 'NR=='$lines' { print $1, $2=$2-2, $3, $4=$4+1, $5, $6=$6-1, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 9 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4=$4+2, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 10 ] || [ "$lines" == 11 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6=$6+2, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi

elif [ "$ru" == "BA_cbond" ]; then
if [ "$lines" == 4 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 5 ] || [ "$lines" == 6 ] || [ "$lines" == 7 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 9 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 11 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
fi

elif [ "$ru" == "BA_nbond" ]; then
if [ "$lines" == 4 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 5 ] || [ "$lines" == 6 ] || [ "$lines" == 7 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6='$anomer', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
elif [ "$lines" == 9 ]; then
awk 'NR=='$lines' { print $1, $2, $3, $4='$anomer', $5, $6='$hemiacetal_O', $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
else
awk 'NR=='$lines' { print $1, $2, $3, $4, $5, $6, $7 }' no_int_coord_base_labels_fixed_v2.txt >> no_int_coord_base_labels_fixed_v3.txt
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

mv internal_coordinates_base_sinH_before.gzmat internal_coordinates_values_base.txt

# Add new lines with new values of distances, angles and dihedrals to file with internal coordinates (this is where I add the numbers of internal coordinates to keep in all molecules)*****IMPORTANT!

# This one allways remains the same

if [ "$type_anomer" == "alpha" ];

then

if [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ];

then

# Insert a1 and d1 for purines
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 98.26" >> first_int_coord_base.txt
echo "d1= 99.56" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for purines
sed -i '2 i a2= 132.67' internal_coordinates_values_base.txt
sed -i '3 i d2= 216.01' internal_coordinates_values_base.txt
sed -i '6 i d3= 185.94' internal_coordinates_values_base.txt

elif [ "$ru" == "cytosine" ] || [ "$ru" == "thymine" ] || [ "$ru" == "uracil" ];

then

# Insert a1 and d1 for pyrimidines
echo "r1= 1.52" > first_int_coord_base.txt
echo "a1= 111.67" >> first_int_coord_base.txt
echo "d1=  97.87" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 109.63' internal_coordinates_values_base.txt
sed -i '3 i d2=  45.02' internal_coordinates_values_base.txt
sed -i '6 i d3= 354.43' internal_coordinates_values_base.txt

elif [ "$ru" == "CA" ];

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

elif [ "$ru" == "melamine" ];

then

# Insert a, d for non canonical base
echo "r1= 2.4765" > first_int_coord_base.txt
echo "a1= 101.43" >> first_int_coord_base.txt # a1=113.40
echo "d1= 117.16" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2=  84.98' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2=  47.71' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 174.76' internal_coordinates_values_base.txt # d3=0.40

sed -i '11s/.*/a5= 174.13/' internal_coordinates_values_base.txt
sed -i '12s/.*/d5= 118.53/' internal_coordinates_values_base.txt

sed -i '20s/.*/a8= 126.13/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 150.03/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9= 106.96/' internal_coordinates_values_base.txt

sed -i '32s/.*/a12=  32.67/' internal_coordinates_values_base.txt
sed -i '33s/.*/d12= 239.28/' internal_coordinates_values_base.txt

sed -i '35s/.*/a13= 148.44/' internal_coordinates_values_base.txt
sed -i '36s/.*/d13= 218.59/' internal_coordinates_values_base.txt

sed -i '38s/.*/a14= 118.14/' internal_coordinates_values_base.txt
sed -i '39s/.*/d14= 204.08/' internal_coordinates_values_base.txt

elif [ "$ru" == "BA_cbond" ]; then

# Insert a, d for BA_cbond
echo "r1= 1.9822" > first_int_coord_base.txt
echo "a1= 142.56" >> first_int_coord_base.txt # a1=113.40
echo "d1=  68.89" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_cbond 
sed -i '2 i a2= 105.94' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 138.68' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3=  83.17' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for BA_cbond
sed -i '8s/.*/a4=  84.85/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 221.94/' internal_coordinates_values_base.txt

sed -i '12s/.*/d5=  73.46/' internal_coordinates_values_base.txt

sed -i '15s/.*/d6= 253.46/' internal_coordinates_values_base.txt

sed -i '18s/.*/d7=  73.44/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9=  87.86/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9=  18.19/' internal_coordinates_values_base.txt

sed -i '30s/.*/d11= 263.16/' internal_coordinates_values_base.txt

elif [ "$ru" == "BA_nbond" ]; then

# Insert a, d for BA_nbond
echo "r1= 1.9919" > first_int_coord_base.txt
echo "a1= 136.98" >> first_int_coord_base.txt # a1=113.40
echo "d1=  68.68" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_nbond 
sed -i '2 i a2=  76.91' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 246.67' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3=  257.47' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for BA_nbond
sed -i '8s/.*/a4= 100.98/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 164.54/' internal_coordinates_values_base.txt

sed -i '12s/.*/d5= 255.59/' internal_coordinates_values_base.txt

sed -i '15s/.*/d6= 75.60/' internal_coordinates_values_base.txt

sed -i '18s/.*/d7= 255.60/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9=  85.96/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9=   8.58/' internal_coordinates_values_base.txt



fi



else # This is where the values of the r, a, d in internal coord of bases for beta sugars are changed



if [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ];

then

# Insert a1 and d1 for purines
echo "a1= 117.78" >> first_int_coord_base.txt
echo "d1= 223.12" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for purines
sed -i '2 i a2= 131.20' internal_coordinates_values_base.txt
sed -i '3 i d2= 334.71' internal_coordinates_values_base.txt
sed -i '6 i d3= 177.00' internal_coordinates_values_base.txt

elif [ "$ru" == "cytosine" ] || [ "$ru" == "thymine" ] || [ "$ru" == "uracil" ];

then

# Insert a1 and d1 for pyrimidines
echo "a1= 108.07" >> first_int_coord_base.txt
echo "d1= 201.62" >> first_int_coord_base.txt

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 116.84' internal_coordinates_values_base.txt
sed -i '3 i d2= 170.38' internal_coordinates_values_base.txt
sed -i '6 i d3= 359.48' internal_coordinates_values_base.txt

elif [ "$ru" == "CA" ];

then

# Insert a, d for non canonical base
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

elif [ "$ru" == "melamine" ];

then

# Insert a, d for non canonical base
echo "r1= 2.4635" > first_int_coord_base.txt
echo "a1= 137.10" >> first_int_coord_base.txt # a1=113.40
echo "d1= 190.69" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for pyrimidines 
sed -i '2 i a2= 84.69' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 188.76' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 173.72' internal_coordinates_values_base.txt # d3=0.40

sed -i '11s/.*/a5= 93.06/' internal_coordinates_values_base.txt
sed -i '12s/.*/d5= 359.75/' internal_coordinates_values_base.txt

sed -i '20s/.*/a8= 126.12/' internal_coordinates_values_base.txt
sed -i '21s/.*/d8= 0.20/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9= 150.03/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9= 182.63/' internal_coordinates_values_base.txt

sed -i '27s/.*/d10= 347.10/' internal_coordinates_values_base.txt

sed -i '30s/.*/d11= 188.57/' internal_coordinates_values_base.txt

sed -i '32s/.*/a12= 33.18/' internal_coordinates_values_base.txt
sed -i '33s/.*/d12= 21.91/' internal_coordinates_values_base.txt

sed -i '34s/.*/r13= 1.3405/' internal_coordinates_values_base.txt
sed -i '35s/.*/a13= 148.56/' internal_coordinates_values_base.txt
sed -i '36s/.*/d13= 358.16/' internal_coordinates_values_base.txt

sed -i '38s/.*/a14= 118.10/' internal_coordinates_values_base.txt
sed -i '39s/.*/d14= 205.80/' internal_coordinates_values_base.txt

elif [ "$ru" == "BA_cbond" ]; then

# Insert a, d for BA_cbond
echo "r1= 2.5763" > first_int_coord_base.txt
echo "a1=  81.20" >> first_int_coord_base.txt # a1=113.40
echo "d1= 197.51" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_cbond 
sed -i '2 i a2= 144.64' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 171.62' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 358.82' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for BA_cbond
sed -i '8s/.*/a4=  59.68/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 170.26/' internal_coordinates_values_base.txt

sed -i '12s/.*/d5= 359.21/' internal_coordinates_values_base.txt

sed -i '15s/.*/d6= 179.20/' internal_coordinates_values_base.txt

sed -i '18s/.*/d7= 359.21/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9=  95.51/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9= 350.59/' internal_coordinates_values_base.txt

sed -i '30s/.*/d11= 178.82/' internal_coordinates_values_base.txt

elif [ "$ru" == "BA_nbond" ]; then

# Insert a, d for BA_nbond
echo "r1= 2.5444" > first_int_coord_base.txt
echo "a1= 137.41" >> first_int_coord_base.txt # a1=113.40
echo "d1= 212.39" >> first_int_coord_base.txt # d1=89.74

# Insert a2 d2 and d3 for BA_nbond 
sed -i '2 i a2=  31.13' internal_coordinates_values_base.txt # a2=116.51
sed -i '3 i d2= 343.58' internal_coordinates_values_base.txt # d2=43.05
sed -i '6 i d3= 179.76' internal_coordinates_values_base.txt # d3=0.40

# Insert other dihedral angles for BA_nbond
sed -i '8s/.*/a4= 116.09/' internal_coordinates_values_base.txt
sed -i '9s/.*/d4= 343.85/' internal_coordinates_values_base.txt

sed -i '12s/.*/d5= 179.86/' internal_coordinates_values_base.txt

sed -i '15s/.*/d6= 359.86/' internal_coordinates_values_base.txt

sed -i '18s/.*/d7= 179.86/' internal_coordinates_values_base.txt

sed -i '23s/.*/a9=  88.71/' internal_coordinates_values_base.txt
sed -i '24s/.*/d9= 163.79/' internal_coordinates_values_base.txt





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
if [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ] || [ "$ru" == "cytosine" ]; then
awk 'NR=='$line' {$0="'$value'= 4.26"} 1' internal_coordinates_base_sinH.gzmat > internal_coordinates_base_sinH_v2.gzmat
elif [ "$ru" == "thymine" ] || [ "$ru" == "uracil" ]; then
awk 'NR=='$line' {$0="'$value'= 181.26"} 1' internal_coordinates_base_sinH.gzmat > internal_coordinates_base_sinH_v2.gzmat
fi


################### Now paste all together to create final z matriz file for the nucleosides ######################
if [ "$ru" == "adenine" ] || [ "$ru" == "guanine" ] || [ "$ru" == "cytosine" ] || [ "$ru" == "thymine" ] || [ "$ru" == "uracil" ]; then
cat no_int_coord_sugar_sinOH.gzmat no_int_coord_base_sinH.gzmat internal_coordinates_values_sugar_sinOH.txt internal_coordinates_base_sinH_v2.gzmat > nucleoside.gzmat
elif [ "$ru" == "CA" ] || [ "$ru" == "melamine" ] || [ "$ru" == "BA_cbond" ] || [ "$ru" == "BA_nbond" ]; then
cat no_int_coord_sugar_sinOH.gzmat no_int_coord_base_sinH.gzmat internal_coordinates_values_sugar_sinOH.txt internal_coordinates_base_sinH.gzmat > nucleoside.gzmat
fi

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


