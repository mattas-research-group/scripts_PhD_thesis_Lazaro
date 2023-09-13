 done#!/bin/bash

# Edit solvated.arc
sed -e '/SUMMARY/,/FINAL/d' solvated.arc > solvated2.mop
sed -i 's/pm7 gnorm=1.0 ef let xyz/pm7 1scf int/g' solvated2.mop 
dos2unix solvated2.mop
grep -v "^$" solvated2.mop > solvated.arc


# Count how many lines are in between pm7
awk '!/pm7/{count++}/pm7/{print count; count = 0}' solvated.arc > linessss.txt
# Eliminate first line
sed '/^$/d' linessss.txt > linessss2.txt
# Eliminate all lines except first in linesssss.txt
sed '2,$d' linessss2.txt > linessss3.txt


# Edit output.q3:
# Eliminate lines before the ones I need
sed -e '/Q3/,/Molecular functions ordered by de(i):/d' output.q3 > temp1.dat
# Eliminate headings of temp1.dat
sed '1,5d' temp1.dat > temp2.dat
# Print 7 column
awk '{ print $7 }' temp2.dat > temporal2.dat
awk '{ print $1*100 }' temporal2.dat > temporal3.dat
# Accumulate populations
sed 'a+p' temporal3.dat | dc -e0 - > porcientos.dat
# Round populations in %
awk 'BEGIN { OFS=FS="\t" } { sub("\\..*", "", $1); print }' porcientos.dat > porcientos_round.dat
# Count populations under 80:
awk '$1 >= 0 && $1 <= 50 { count++ } END { print count ? count : 0 }' porcientos_round.dat > count.dat
# Add 1 to count.dat
a=`cat count.dat`
export b=$((1+$a))

# Select cells with the population
awk '{print $1}' temp2.dat > temp3.dat
head -$b temp3.dat | tail -$b > temp4.dat
sort -n temp4.dat > list.dat

# Select cells
w=`cat linessss3.txt`
cat list.dat | while read n; do 
awk "/CELL:[ ]*$n$/ {print NR}" solvated.arc >> n_lines.dat
done 
while read m; do 
h=$(($m-2))
echo $h >> lineas.dat
i=$(($h+$w))
sed -n "$h,$i p"  solvated.arc >> solvated.mop
done <n_lines.dat
sed -i 's/pm7 1scf int/pm7 gnorm=1.0 ef let xyz/g' solvated.mop
rm solvated2.mop temp1.dat temp2.dat temp3.dat n_lines.dat lineas.dat temporal2.dat temporal3.dat porcientos.dat count.dat linessss.txt linessss2.txt 





