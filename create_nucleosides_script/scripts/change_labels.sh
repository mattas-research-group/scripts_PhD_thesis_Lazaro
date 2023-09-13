#!/bin/bash

# This script swich the labels between two hydrogen atoms in nucleoside to apply the script to create nucleotide 

# Insert the names of the files for the TC and the RUs
read -p "insert original label of hydrogen that will react:" old_label
read -p "insert new label for hydrogen that will react:" new_label

# Substitute the labels of the hydrogen to be replaced to create the nucleotides
line_n_old_label=$(($old_label+5))
line_n_new_label=$(($new_label+5))

# Separate file in internal coordinates values files and z matrix file:
sed -n '/Variables:/,$p' dA_fura_2endo_alpha.gzmat > internal_coordinates_values.txt
sed '1,/Variables:/!d' dA_fura_2endo_alpha.gzmat > z_matrix.txt
# Delete last line of this file 
sed -i '$d' z_matrix.txt

# Swap the two lines in the z matrix file (e.g. put old line in new and biseversa)
printf '%s\n' "$line_n_old_label"m"$line_n_new_label" "$line_n_new_label"-m"$line_n_old_label"- w q | ed -s z_matrix.txt

# Replace rlabels, alabels and dlabels in Z matrix of nucleoside
sed -i 's/r'$old_label'/r'$new_label'/g' z_matrix.txt
sed -i 's/a'$old_label'/a'$new_label'/g' z_matrix.txt 
sed -i 's/d'$old_label'/d'$new_label'/g' z_matrix.txt

# Replace fist instance of new labels by old labels
sed -i ':a;N;$!ba;s/r'$new_label'/r'$old_label'/1' z_matrix.txt
sed -i ':a;N;$!ba;s/a'$new_label'/a'$old_label'/1' z_matrix.txt
sed -i ':a;N;$!ba;s/d'$new_label'/d'$old_label'/1' z_matrix.txt

# Obtain number of lines in internal coordinates values files for react hydrogen
number_line_r_old_label=$(($old_label*3-7))
number_line_a_old_label=$(($old_label*3-6))
number_line_d_old_label=$(($old_label*3-5))
echo $number_line_d_old_label

# Obtain number of lines in internal coordinates values files for non-react hydrogen
number_line_r_new_label=$(($new_label*3-7))
number_line_a_new_label=$(($new_label*3-6))
number_line_d_new_label=$(($new_label*3-5))
echo $number_line_d_new_label

# Swap lines that contain new and old rlabel, alabel and dlabel
printf '%s\n' "$number_line_r_old_label"m"$number_line_r_new_label" "$number_line_r_new_label"-m"$number_line_r_old_label"- w q | ed -s internal_coordinates_values.txt 
printf '%s\n' "$number_line_a_old_label"m"$number_line_a_new_label" "$number_line_a_new_label"-m"$number_line_a_old_label"- w q | ed -s internal_coordinates_values.txt 
printf '%s\n' "$number_line_d_old_label"m"$number_line_d_new_label" "$number_line_d_new_label"-m"$number_line_d_old_label"- w q | ed -s internal_coordinates_values.txt 

# Replace rlabels, alabels and dlabels in internal coord values file of nucleoside
sed -i 's/r'$old_label'=/r'$new_label'=/g' internal_coordinates_values.txt 
sed -i 's/a'$old_label'=/a'$new_label'=/g' internal_coordinates_values.txt
sed -i 's/d'$old_label'=/d'$new_label'=/g' internal_coordinates_values.txt

sed -i ':a;N;$!ba;s/r'$new_label'=/r'$old_label'=/1' internal_coordinates_values.txt
sed -i ':a;N;$!ba;s/a'$new_label'=/a'$old_label'=/1' internal_coordinates_values.txt
sed -i ':a;N;$!ba;s/d'$new_label'=/d'$old_label'=/1' internal_coordinates_values.txt

# Put z matrix and internal coord values together
cat z_matrix.txt internal_coordinates_values.txt > new.gzmat

# Convert gzmat to mol format using obabel
obabel new.gzmat -O new.mol



