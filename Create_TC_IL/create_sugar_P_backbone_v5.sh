#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a IL_names=("phosphate" "arsenate")


function create_TC_IL_backbone {

# Note1: no_int_coord_base.gzmat and no_int_coord_sugar.gzmat are the files with only the z matrix
# Note2: this initial test was with beta 2deoxyribofuranose and guanine

# Insert the names of the files for the TC and the RUs
#read -p "insert environment of nucleoside (e.g. vacuum or water)" environment
#read -p "insert name of log file for the Ionized linker (IL) without extension (e.g. phosphate or arsenate):" IL
#read -p "insert type of anomer in case that it is (e.g. beta, or alpha, or none if there is not anomer)" type_anomer
#read -p "insert type of trifunctional connector (TC) in case that it is (e.g. 2deoxy_ribofuranose, or ribofuranose, or 2deoxy_ribopyranose, or ribopyranose, or threose, or glycerol, or peptide)" type_sugar
#read -p "insert conformation of the sugar ring (e.g. 2endo_alpha, or 2endo_beta, or D1C4_alpha, or D4C1_beta, or none if it is glycerol, glyceric acid or peptide)" sugar_ring

# Declare corresponding variables
local environment=$1
local type_sugar=$2
local sugar_ring=$3
local IL=$4

# Give name to trifunctional connector file
if [ "$type_sugar" == "glycerol" ] || [ "$type_sugar" == "glyceric_acid" ]; then
tc=$type_sugar
else
tc="$type_sugar"_"$sugar_ring"
fi

# Obtain type of anomer
if [ "$type_sugar" != "glycerol" ] && [ "$type_sugar" != "glyceric_acid" ]; then
cut -d "_" -f2 <<< "$sugar_ring" > type_anomer.txt
type_anomer=`cat type_anomer.txt`
fi

# Convert tc to gzmat
if [ "$sugar_ring" == "D1C4_beta" ] || [ "$sugar_ring" == "D4C1_beta" ]; then
obabel $tc.mol -O $tc.gzmat
obabel $tc.mol -O $tc.pdb
else
obabel $tc.log -O $tc.gzmat
obabel $tc.log -O $tc.pdb
fi

# Convert IL to gzmat
obabel $IL.log -O $IL.gzmat

################ Insert options to authomatically change labels of atoms to delete and hemiacetal O and anomeric carbon for sugar ########################

if [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "2deoxy_ribofuranose" ]; then 
if [ "$sugar_ring" == "2endo_beta" ] || [ "$sugar_ring" == "2endo_alpha" ]; then
anomer=1
hemiacetal_O=2
hydrogen=9
elif [ "$sugar_ring" == "3endo_beta" ] || [ "$sugar_ring" == "3endo_alpha" ]; then
anomer=1
hemiacetal_O=2
hydrogen=16
fi

elif [ "$type_sugar" == "2deoxy_ribopyranose" ]; then
if [ "$sugar_ring" == "D1C4_alpha" ] || [ "$sugar_ring" == "D4C1_alpha" ]; then
anomer=14
hemiacetal_O=4
hydrogen=15
elif [ "$sugar_ring" == "D1C4_beta" ] || [ "$sugar_ring" == "D4C1_beta" ]; then
anomer=13
hemiacetal_O=4
hydrogen=14
fi

elif [ "$type_sugar" == "ribopyranose" ]; then 
if [ "$sugar_ring" == "D1C4_alpha" ] || [ "$sugar_ring" == "D4C1_alpha" ]; then
anomer=16
hemiacetal_O=4
hydrogen=17
elif [ "$sugar_ring" == "D1C4_beta" ] || [ "$sugar_ring" == "D4C1_beta" ]; then
anomer=13
hemiacetal_O=4
hydrogen=14
fi

elif [ "$type_sugar" == "threose" ]; then
anomer=13
hemiacetal_O=3
hydrogen=14

elif [ "$type_sugar" == "glycerol" ]; then

anomer=5
hemiacetal_O=4
hydrogen=11

elif [ "$type_sugar" == "glyceric_acid" ]; then
anomer=4
hemiacetal_O=3
hydrogen=8



fi

echo "the anomer atom is $anomer and the hydrogen to delete from TC is $hydrogen"

############################################################### Create if then statement to work either in peptide or other TC ########################################

########################
## WORK ON SUGAR FILE ##
########################

if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ]; then
atom_before_hemiacetal=$(($hemiacetal_O+1))
atom_before_before_hemiacetal=$(($hemiacetal_O+2))

elif [ "$type_sugar" == "threose" ]; then
atom_before_hemiacetal=$(($hemiacetal_O-2))
atom_before_before_hemiacetal=$(($atom_before_hemiacetal+1))

elif [ "$type_sugar" == "glycerol" ]; then
atom_before_hemiacetal=$(($hemiacetal_O-2))

elif [ "$type_sugar" == "glyceric_acid" ]; then
atom_before_hemiacetal=$(($hemiacetal_O-2))

fi

# Obtain the lines number for these atoms
line_hydrogen=$(($hydrogen+6))

# Delete lines in matrix that correspond to the O and H in anomeric position
sed ''$line_hydrogen','$line_hydrogen'd' $tc.gzmat > sugar_sin_H.gzmat

#### Delete in internal coordinates values lines with values of the distances, angles and dihedrals for the atoms deleted in previous step ###
sed -i '/r'$hydrogen'/d' sugar_sin_H.gzmat
sed -i '/a'$hydrogen'/d' sugar_sin_H.gzmat
sed -i '/d'$hydrogen'/d' sugar_sin_H.gzmat

########### Calculate total number of atoms in sugar ############

# Delete first headings of this file
sed '1,6d' $tc.gzmat > sinheadings_sugar.gzmat

# Delete values the internal coordinates
sed -i '1,/Variables:/!d' sinheadings_sugar.gzmat
sed -i '$d' sinheadings_sugar.gzmat

# Count number of lines in this file
< sinheadings_sugar.gzmat wc -l > number_atoms_sugar.txt

atoms_sugar=`cat number_atoms_sugar.txt`
atoms_sugar_ion=$(($atoms_sugar-1))
lines_plus_atoms_in_sugar_ion=$(($atoms_sugar_ion+6))
starting_line_int_coord=$(($atoms_sugar_ion+7))

#################################################################

# Copy lines with values of internal coordinates to another file
sed -n ''$starting_line_int_coord',$p' sugar_sin_H.gzmat > internal_coordinates_values_sugar_sinH.txt 
# Eliminate last empty line
sed -i '$ d' internal_coordinates_values_sugar_sinH.txt      

# Delete values the internal coordinates
sed -i '1,/Variables:/!d' sugar_sin_H.gzmat

# Delete last line of this file 
sed -i '$d' sugar_sin_H.gzmat

### Change r#, a# and d# and atoms labels in file of sugar from line 4 by new numbers

#line_hydrogen=$(($hydrogen+6))

for lines in $(seq 7 1 $lines_plus_atoms_in_sugar_ion)

do

if [[ "$lines" -ge "$line_hydrogen" ]]; then
# Calculate new r# a# and d#
label_to_replace=$(($lines-5))
new_label=$(($label_to_replace-1))
sed -i ''$lines's/r'$label_to_replace'/r'$new_label'/' sugar_sin_H.gzmat
sed -i ''$lines's/a'$label_to_replace'/a'$new_label'/' sugar_sin_H.gzmat
sed -i ''$lines's/d'$label_to_replace'/d'$new_label'/' sugar_sin_H.gzmat
fi

done

# Copy the initial instructions in final file
sed -n '1,6p' sugar_sin_H.gzmat > instructions_sugar_sin_H.gzmat

# Eliminate instructions from this file
sed '1,6d' sugar_sin_H.gzmat > sugar_sin_H_sin_instructions.gzmat

# Correct atoms labels for the case of missing a hydrogen
for column in $(seq 1 1 7); do
awk '{ print $'$column' }' sugar_sin_H_sin_instructions.gzmat > column_"$column"_from_sugar_sin_H.txt

if [[ "$column" -eq 2 ]] || [[ "$column" -eq 4 ]] || [[ "$column" -eq 6 ]]; then

cat column_"$column"_from_sugar_sin_H.txt | while read atom; 

do

if [[ "$atom" -lt "$hydrogen" ]]; then
echo $atom >> column_"$column"_from_sugar_sin_H2.txt

elif [[ "$atom" -gt "$hydrogen" ]]; then
new_atom=$(($atom-1))
echo $new_atom >> column_"$column"_from_sugar_sin_H2.txt

fi

done

rm column_"$column"_from_sugar_sin_H.txt
mv column_"$column"_from_sugar_sin_H2.txt column_"$column"_from_sugar_sin_H.txt

fi
done

# Put all the columns together
paste column_1_from_sugar_sin_H.txt column_2_from_sugar_sin_H.txt column_3_from_sugar_sin_H.txt column_4_from_sugar_sin_H.txt column_5_from_sugar_sin_H.txt column_6_from_sugar_sin_H.txt column_7_from_sugar_sin_H.txt > temp.txt

# Reduce spaces between columns in temp.txt
sed -i $'s/\t/  /g' temp.txt

# Add instructions
cat instructions_sugar_sin_H.gzmat temp.txt > no_int_coord_sugar_sinH.gzmat ### THIS IS THE FINAL FILE FOR THE Z MATRIX OF SUGAR WITHOUT THE INT COORD VALUES


############# Work on correcting the rlabel, alabel and dlabel in the file without internal coordinates values ################

for labels in $(seq 2 1 $atoms_sugar_ion)

do

echo "r"$labels"=" >> labels_sugar_fix.txt
echo "a"$labels"=" >> labels_sugar_fix.txt
echo "d"$labels"=" >> labels_sugar_fix.txt

done

# Eliminate lines from the previous file that correspond to a2, d2 and d3
sed -i '2,3d' labels_sugar_fix.txt
sed -i '4d' labels_sugar_fix.txt

# Remove heading
sed -i '1,1d' internal_coordinates_values_sugar_sinH.txt

# Eliminate everything before the = sign in old int coord file from the base and eliminate first line
cut -f2 -d"=" internal_coordinates_values_sugar_sinH.txt > internal_coordinates_values_sugar_bef_nucleo.txt

# Paste new labels and values for internal coordinates
paste labels_sugar_fix.txt internal_coordinates_values_sugar_bef_nucleo.txt > new_internal_coordinates_sugar_bef_nucleo.txt

# Reduce spaces between two columns
sed -i $'s/\t//g' new_internal_coordinates_sugar_bef_nucleo.txt
tr -s " " < new_internal_coordinates_sugar_bef_nucleo.txt > internal_coordinates_sugar_sinH_before.gzmat

# Delete previous file
rm internal_coordinates_values_sugar_sinH.txt

# Add heading to internal coordinates values
sed -i '1 i Variables:' internal_coordinates_sugar_sinH_before.gzmat

cp internal_coordinates_sugar_sinH_before.gzmat internal_coordinates_values_sugar_sinH.gzmat ### THIS IS THE FINAL FILE FOR INTERNAL COORDINATES VALUES OF THE SUGAR

rm sinheadings_sugar.gzmat number_atoms_sugar.txt 

#######################
## WORK ON ILs FILE  ##
#######################

# Insert labels of atoms to delete for IL 
oxygen=4
hydrogen=5

# Delete headings from base file
sed -i '1,6d' $IL.gzmat

### Calculate number of atoms in base #######

# Delete values the internal coordinates
sed '1,/Variables:/!d' $IL.gzmat > no_int_coord_IL.txt

# Delete last line of this file 
sed -i '$d' no_int_coord_IL.txt

# Count number of lines in this file and determine number of atoms in IL
< no_int_coord_IL.txt wc -l > number_atoms_IL_new.txt

atoms_IL=`cat number_atoms_IL_new.txt`
echo $atoms_IL
ions_atoms_IL=$(($atoms_IL-2))
echo $ions_atoms_IL
first_atom_IL=$(($atoms_sugar_ion+1))
second_atom_IL=$(($atoms_sugar_ion+2))
third_atom_IL=$(($atoms_sugar_ion+3))
fourth_atom_IL=$(($atoms_sugar_ion+4))
five_atom_IL=$(($atoms_sugar_ion+5))
#seven_atom_base=$(($atoms_sugar_ion+7))
#nine_atom_base=$(($atoms_sugar_ion+9))

# Delete oxygen and hydrogen from IL Z-matrix 
sed -i ''$oxygen','$hydrogen'd' no_int_coord_IL.txt

### Change r#, a# and d# and atoms labels in file of IL from line 4 by new numbers

for lines in $(seq 4 1 $ions_atoms_IL)

do

# Calculate new r# a# and d#
new_lines=$(($lines+2))
new_labels_number_IL=$(($atoms_sugar_ion+$lines))
sed -i ''$lines's/r'$new_lines'/r'$new_labels_number_IL'/' no_int_coord_IL.txt
sed -i ''$lines's/a'$new_lines'/a'$new_labels_number_IL'/' no_int_coord_IL.txt
sed -i ''$lines's/d'$new_lines'/d'$new_labels_number_IL'/' no_int_coord_IL.txt


# Calculate new atoms labels for the 
awk 'NR=='$lines' { print $1,  $2=$2+'$atoms_sugar_ion',  $3,  $4=$4+'$atoms_sugar_ion',  $5,  $6=$6+'$atoms_sugar_ion',  $7 }' no_int_coord_IL.txt >> no_int_coord_IL_labels_fixed.txt 

done

## Add the three first lines to the z matriz of the base and fix them accordingly
# Note: to do this is important copy first column of this file and three first rows in independent files
awk '{ print $1 }' no_int_coord_IL.txt > atoms_symbols_no_int_coord_IL.txt

# Print first, second and third symbols
sed -n 1,1p atoms_symbols_no_int_coord_IL.txt > first_symbol.txt
first_symbol=`cat first_symbol.txt`
sed -n 2,2p atoms_symbols_no_int_coord_IL.txt > second_symbol.txt
second_symbol=`cat second_symbol.txt`
sed -n 3,3p atoms_symbols_no_int_coord_IL.txt > third_symbol.txt
third_symbol=`cat third_symbol.txt`

# Replace first line with internal coordinates as it should be in file
sed -i '1 s/.*/\'$first_symbol'  '$anomer'  r'$first_atom_IL'  '$hemiacetal_O'  a'$first_atom_IL'  '$atom_before_hemiacetal'  d'$first_atom_IL'/' no_int_coord_IL.txt

# Replace second line with another specific coordinates
sed -i '2 s/.*/\'$second_symbol'  '$first_atom_IL'  r'$second_atom_IL'  '$anomer'  a'$second_atom_IL'  '$hemiacetal_O'  d'$second_atom_IL'/' no_int_coord_IL.txt

# Replace third line with another specific coordinates
sed -i '3 s/.*/\'$third_symbol'  '$second_atom_IL'  r'$third_atom_IL'  '$first_atom_IL'  a'$third_atom_IL'  '$anomer'  d'$third_atom_IL'/' no_int_coord_IL.txt

# Replace the rest of the 

# Delete lines from line 3 in no_int_coord_base.txt
sed -i '4,$d' no_int_coord_IL.txt

# Paste all together 
cat no_int_coord_IL.txt no_int_coord_IL_labels_fixed.txt > no_int_coord_IL_labels_fixed_v2.txt

# Replace rest of lines
for lines in $(seq 1 1 $ions_atoms_IL); do

#if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ] || [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ] || [ "$type_sugar" == "threose" ] || [ "$type_sugar" == "glycerol" ] || [ "$type_sugar" == "glyceric_acid" ]; then

if [ "$lines" == 4 ] || [ "$lines" == 5 ]; then
awk 'NR=='$lines' { print $1,  $2,  $3,  $4='$anomer',  $5,  $6='$hemiacetal_O',  $7 }' no_int_coord_IL_labels_fixed_v2.txt >> no_int_coord_IL_labels_fixed_v3.txt
else 
awk 'NR=='$lines' { print $1,  $2,  $3,  $4,  $5,  $6,  $7 }' no_int_coord_IL_labels_fixed_v2.txt >> no_int_coord_IL_labels_fixed_v3.txt
fi

#elif

#fi

done

# Rename this file as final
mv no_int_coord_IL_labels_fixed_v3.txt no_int_coord_IL_sinOH.gzmat  #### THIS IS THE FINAL FILE FOR THE Z MATRIZ OF THE IL FIXED APPROPIATELY

################ Work in internal coordinates values of IL ############################

# Copy values of internal coordinates in different file
lines2=$(($ions_atoms_IL+1))

# Copy lines with values of internal coordinates to another file
sed -n '/Variables:/,$p' $IL.gzmat > internal_coordinates_values_IL.txt

# Delete first line of this file
sed -i '1d' internal_coordinates_values_IL.txt

# Delete coordinates values for oxygen and hydrogen that was deleted from IL
sed -i '/r'$oxygen'/d' internal_coordinates_values_IL.txt
sed -i '/a'$oxygen'/d' internal_coordinates_values_IL.txt
sed -i '/d'$oxygen'/d' internal_coordinates_values_IL.txt
sed -i '/r'$hydrogen'/d' internal_coordinates_values_IL.txt
sed -i '/a'$hydrogen'/d' internal_coordinates_values_IL.txt
sed -i '/d'$hydrogen'/d' internal_coordinates_values_IL.txt

# Change labels of each atoms in internal coordinates file of RU before attaching to TC

for labels in $(seq 1 1 $ions_atoms_IL)

do

echo "r$labels=" >> labels_IL_fix.txt
echo "a$labels=" >> labels_IL_fix.txt
echo "d$labels=" >> labels_IL_fix.txt

done

# Insert line in line 2, 3 and 6
# sed -i 
#sed -i '2i a2=' internal_coordinates_values_IL.txt
#sed -i '3i d2=' internal_coordinates_values_IL.txt
#sed -i '6i d3=' internal_coordinates_values_IL.txt

# Eliminate everything before the = sign in old int coord file from the base
cut -f2 -d"=" internal_coordinates_values_IL.txt > internal_coordinates_values_IL_bef_nucleo.txt

# Paste new labels and values for internal coordinates
paste labels_IL_fix.txt internal_coordinates_values_IL_bef_nucleo.txt > new_internal_coordinates_IL_bef_nucleo.txt

# Reduce spaces between two columns
sed -i $'s/\t//g' new_internal_coordinates_IL_bef_nucleo.txt
tr -s " " < new_internal_coordinates_IL_bef_nucleo.txt > internal_coordinates_IL_sinOH_before.gzmat

# Delete previous file
rm internal_coordinates_values_IL.txt

cp internal_coordinates_IL_sinOH_before.gzmat internal_coordinates_values_IL.txt

# Add new lines with new values of distances, angles and dihedrals to file with internal coordinates (this is where I add the numbers of internal coordinates to keep in all molecules)*****IMPORTANT!

#if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ]; then
# Insert a, d for non canonical base
echo "r1=   1.6" > first_int_coord_IL.txt
echo "a1= 119.01" >> first_int_coord_IL.txt # a1=113.40
echo "d1= 179.97" >> first_int_coord_IL.txt # d1=89.74

# Insert a2 d2 and d3  
sed -i '2 i a2= 100.68' internal_coordinates_values_IL.txt # a2=116.51
sed -i '3 i d2= 172.83' internal_coordinates_values_IL.txt # d2=43.05
sed -i '6 i d3= 269.96' internal_coordinates_values_IL.txt # d3=0.40

sed -i '8s/.*/a4= 106.20/' internal_coordinates_values_IL.txt
sed -i '9s/.*/d4= 283.85/' internal_coordinates_values_IL.txt

sed -i '11s/.*/a5= 107.45/' internal_coordinates_values_IL.txt
sed -i '12s/.*/d5=  61.25/' internal_coordinates_values_IL.txt



# Put both files together
cat first_int_coord_IL.txt internal_coordinates_values_IL.txt > internal_coordinates_IL.txt

# Change labels of each atoms in internal coordinates file of IL

for labels in $(seq 1 1 $ions_atoms_IL)

do

new_label=$(($labels+$atoms_sugar_ion))

echo "r$new_label=" >> newlabels_IL.txt
echo "a$new_label=" >> newlabels_IL.txt
echo "d$new_label=" >> newlabels_IL.txt

done

# Eliminate everything before the = sign in old int coord file from the base
cut -f2 -d"=" internal_coordinates_IL.txt > new_internal_coordinates_IL.txt

# Paste new labels and values for internal coordinates
paste newlabels_IL.txt new_internal_coordinates_IL.txt > new_internal_coordinates_IL2.txt

# Reduce spaces between two columns
sed -i $'s/\t//g' new_internal_coordinates_IL2.txt
tr -s " " < new_internal_coordinates_IL2.txt > internal_coordinates_IL_sin_OH.gzmat #### THIS IS THE FINAL FILE WITH THE INTERNAL COORD VALUES FIXED FOR THE BASE 

# Fix value of dihedral angle between first hydrogen atom of base:
cat no_int_coord_sugar_sinH.gzmat no_int_coord_IL_sinOH.gzmat internal_coordinates_values_sugar_sinH.gzmat internal_coordinates_IL_sin_OH.gzmat > backbone.gzmat

# Adjust dihedral that put all CH2HPO4 group together equal #atoms in sugarphospate + number atoms for dihedral*3 + 1
if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ]; then
#if [ "$sugar_ring" == "2endo_alpha" ] || [ "$sugar_ring" == "3endo_beta" ]; then

line_to_change_d4=$(($atoms_sugar_ion+$ions_atoms_IL+13))
line_to_change_d5=$(($atoms_sugar_ion+$ions_atoms_IL+16))
line_to_change_d11=$(($atoms_sugar_ion+$ions_atoms_IL+34))
echo "line of d4 is #$line_to_change_d4"
echo "line of d5 is #$line_to_change_d5"
echo "line of d11 is #$line_to_change_d11"

#sed -i ''$line_to_change_d4's/.*/d4= 275.90/' backbone.gzmat
#sed -i ''$line_to_change_d5's/.*/d5=  30.89/' backbone.gzmat
#sed -i ''$line_to_change_d11's/.*/d11= 155.91/' backbone.gzmat

#elif [ "$sugar_ring" == "2endo_beta" ] || [ "$sugar_ring" == "3endo_alpha" ]; then

#line_to_change_d4=$(($atoms_sugar_ion+$ions_atoms_IL+13))
#line_to_change_d5=$(($atoms_sugar_ion+$ions_atoms_IL+16))
#line_to_change_d11=$(($atoms_sugar_ion+$ions_atoms_IL+34))
#echo "line of d4 is #$line_to_change_d4"
#echo "line of d5 is #$line_to_change_d5"
#echo "line of d11 is #$line_to_change_d11"

#sed -i ''$line_to_change_d4's/.*/d4= 275.90/' backbone.gzmat
#sed -i ''$line_to_change_d5's/.*/d5=  30.89/' backbone.gzmat
#sed -i ''$line_to_change_d11's/.*/d11= 155.91/' backbone.gzmat

#fi

echo "* $hemiacetal_O $atom_before_hemiacetal *" >> coordinates_to_scan.log
echo "D $anomer $hemiacetal_O $atom_before_hemiacetal 5 S 6 60.000000" >> coordinates_to_scan.log 

#elif [ "$type_sugar" == "ribofuranose" ]; then
#if [ "$sugar_ring" == "2endo_alpha" ] || [ "$sugar_ring" == "2endo_beta" ]; then

#line_to_change_d4=$(($atoms_sugar_ion+$ions_atoms_IL+13))
#line_to_change_d5=$(($atoms_sugar_ion+$ions_atoms_IL+16))
#line_to_change_d11=$(($atoms_sugar_ion+$ions_atoms_IL+34))
#echo "line of d4 is #$line_to_change_d4"
#echo "line of d5 is #$line_to_change_d5"
#echo "line of d11 is #$line_to_change_d11"

#sed -i ''$line_to_change_d4's/.*/d4= 275.90/' backbone.gzmat
#sed -i ''$line_to_change_d5's/.*/d5=  30.89/' backbone.gzmat
#sed -i ''$line_to_change_d11's/.*/d11= 155.91/' backbone.gzmat

#elif [ "$sugar_ring" == "3endo_beta" ] || [ "$sugar_ring" == "3endo_alpha" ]; then

#line_to_change_d4=$(($atoms_sugar_ion+$ions_atoms_IL+13))
#line_to_change_d5=$(($atoms_sugar_ion+$ions_atoms_IL+16))
#line_to_change_d11=$(($atoms_sugar_ion+$ions_atoms_IL+34))
#echo "line of d4 is #$line_to_change_d4"
#echo "line of d5 is #$line_to_change_d5"
#echo "line of d11 is #$line_to_change_d11"

#sed -i ''$line_to_change_d4's/.*/d4= 275.90/' backbone.gzmat
#sed -i ''$line_to_change_d5's/.*/d5=  30.89/' backbone.gzmat
#sed -i ''$line_to_change_d11's/.*/d11= 155.91/' backbone.gzmat

#fi

#echo "* $hemiacetal_O $atom_before_hemiacetal *" >> coordinates_to_scan.log
#echo "D $anomer $hemiacetal_O $atom_before_hemiacetal $atom_before_before_hemiacetal S 6 60.000000" >> coordinates_to_scan.log 

elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ]; then

echo "* $anomer $hemiacetal_O *" >> coordinates_to_scan.log
echo "D $first_atom_IL $anomer $hemiacetal_O 3 S 6 60.000000" >> coordinates_to_scan.log 

elif [ "$type_sugar" == "threose" ]; then

echo "* $anomer $hemiacetal_O *" >> coordinates_to_scan.log
echo "D $first_atom_IL $anomer $hemiacetal_O 4 S 6 60.000000" >> coordinates_to_scan.log 

elif [ "$type_sugar" == "glycerol" ] || [ "$type_sugar" == "glyceric_acid" ]; then

echo "* $anomer $hemiacetal_O *" >> coordinates_to_scan.log
echo "D $first_atom_IL $anomer $hemiacetal_O $atom_before_hemiacetal S 6 60.000000" >> coordinates_to_scan.log 

fi

# Replace in final gzmat charge and multiplicity
sed -i '6s/.*/-1 1/' backbone.gzmat

# Create new name of nucleoside
backbone="$tc"_"$IL"

# Convert this file to gjf
obabel backbone.gzmat -O $backbone.gjf

rm *.txt
rm *.gzmat

# Delete first two lines from final gjf file
sed -i 1,2d $backbone.gjf

# Attach DFT instructions
cat $MD_RUN/instructions_DFT/instructions_b3lyp.txt $backbone.gjf > "$backbone"2.gjf

# Attach coordinates to scan
cat "$backbone"2.gjf coordinates_to_scan.log > "$backbone"_final.gjf

# Substitute keywords for the case of water opt
if [ "$environment" == "water" ]; then
sed -i 's/scf=maxcycle=1600/scf=maxcycle=1600 scrf=(solvent=water)/g' "$backbone"_final.gjf
fi

# Replace chk file by corresponding correct name
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$backbone'_scan.chk/g' "$backbone"_final.gjf

# Delete unccesary gjf and log files
rm "$backbone".gjf "$backbone"2.gjf coordinates_to_scan.log

# Change final name of backbone gjf file
mv "$backbone"_final.gjf "$backbone"_scan.gjf

# Copy script to run calculation
cp $MD_RUN/job1.sh .

# Replace name of file by backbone name
sed -i 's/g16 input.gjf/g16 '$backbone'_scan.gjf/g' job1.sh


}

##########################################################################################################################
##########################################################################################################################

############################################## WORK ON EACH FOLDER #######################################################

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# Open corresponding folder
cd TCs_ILs_backbone/gaussian_scan


for environment in "${environment[@]}"; do

#mkdir $environment

cd $environment

for TC in "${TC_names[@]}"; do

cd $TC

if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

for ring_conf_fura in "${Ring_Conf_furanose[@]}"; do

cd $ring_conf_fura

for IL in "${IL_names[@]}"; do

cd $IL

# Delete existing gjf file
rm *.gjf
rm job1.sh

# Run function to generate the TC-IL backbone
create_TC_IL_backbone $environment $TC $ring_conf_fura $IL

# Print message of status
echo "working in $environment environment, TC $TC with $ring_conf_fura pyranose ring conf and IL $IL"


cd ../
done

cd ../
done

elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then

for ring_conf_pyra in "${Ring_Conf_pyranose[@]}"; do

cd $ring_conf_pyra

for IL in "${IL_names[@]}"; do

cd $IL

# Delete existing gjf file
rm *.gjf
rm job1.sh

# Run function to generate the TC-IL backbone
create_TC_IL_backbone $environment $TC $ring_conf_pyra $IL

# Print message of status
echo "working in $environment environment, TC $TC with $ring_conf_pyra pyranose ring conf and IL $IL"


cd ../
done

cd ../
done

elif [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then

for IL in "${IL_names[@]}"; do

cd $IL

ring_conf="none"

# Delete existing gjf file
rm *.gjf
rm job1.sh

# Run function to generate the TC-IL backbone
create_TC_IL_backbone $environment $TC $ring_conf $IL

# Print message of status
echo "working in $environment environment, TC $TC with none sugar ring conf and IL $IL"





cd ../
done

fi

cd ../
done


cd ../
done
