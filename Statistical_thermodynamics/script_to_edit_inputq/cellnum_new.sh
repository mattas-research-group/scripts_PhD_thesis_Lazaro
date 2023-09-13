#!/bin/bash
grep -e "CELL:" solvated.arc |uniq | awk '{FS=OFS=":"} {print $2}'|sed -e 's/ //g' > cell_list 
