#!/bin/bash

lines_plus_atoms_in_sugar_ion=36
line_hydrogen=15
hydrogen=9
atoms_sugar_ion=30
ions_atoms_IL=5
anomer=1
hemiacetal_O=2

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
