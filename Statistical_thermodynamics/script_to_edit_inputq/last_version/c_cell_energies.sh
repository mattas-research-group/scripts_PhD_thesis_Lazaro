#!/bin/bash
grep -e "TOTAL ENERGY" solvated.arc | awk '{print $4}' > cell_energy
grep -e "CELL:" solvated.arc |uniq | awk '{FS=OFS=":"} {print $2}'|sed -e 's/ //g' > cell_list 
awk 'FNR==NR{a[NR]=$1;next}{$5=a[FNR]}1' cell_list cell_energy > temp1
sed -i -e 's/^/    /' temp1
sed -n '1,2p' input.q > 2lines
cat 2lines temp1 > input.q
rm cell_energy cell_list temp1 2lines 



