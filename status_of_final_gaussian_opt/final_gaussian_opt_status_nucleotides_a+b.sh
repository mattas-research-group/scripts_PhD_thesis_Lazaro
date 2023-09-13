#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a RU_names=("adenine" "guanine" "cytosine" "thymine" "uracil" "TARC_cbond" "TARC_nbond" "BA_cbond" "BA_nbond" "CA" "melamine")
declare -a RU_ac_names=("adenine_ac" "guanine_ac" "cytosine_ac" "thymine_ac" "uracil_ac" "TARC_cbond_ac" "TARC_nbond_ac" "BA_cbond_ac" "BA_nbond_ac" "CA_ac" "melamine_ac")
declare -a IL_names=("phosphate" "arsenate")



function completion_status {

# Assign variables names to entries of function
local environment=$1
local type_sugar=$2
local sugar_ring=$3
local IL=$4
local ru=$5
########################################################## Determine if the calculation converged or not #######################################################

# Print number of line of last instance of Maximum Force
if [ "$sugar_ring" == "none" ]; then
awk '/Maximum Force/ { ln = FNR } END { print ln }' "$type_sugar"_"$ru"_"$IL"_min.log > line_maximumforce_"$type_sugar"_"$ru"_"$IL".txt
# Assign a variable to this file 
line=`cat line_maximumforce_"$type_sugar"_"$ru"_"$IL".txt`
else
awk '/Maximum Force/ { ln = FNR } END { print ln }' "$type_sugar"_"$sugar_ring"_"$ru"_"$IL"_min.log > line_maximumforce_"$type_sugar"_"$sugar_ring"_"$ru"_"$IL".txt
# Assign a variable to this file 
line=`cat line_maximumforce_"$type_sugar"_"$sugar_ring"_"$ru"_"$IL".txt`
fi

# Obtain line that contains convergence message
total_lines=$(($line+6)) 

# Copy line with convergence message in new file
if [ "$sugar_ring" == "none" ]; then
sed -n ''$total_lines'p' "$type_sugar"_"$ru"_"$IL"_min.log > convergence_"$type_sugar"_"$ru".txt
# Assign variable to this file
convergence=`cat convergence_"$type_sugar"_"$ru".txt`
else
sed -n ''$total_lines'p' "$type_sugar"_"$sugar_ring"_"$ru"_"$IL"_min.log > convergence_"$type_sugar"_"$sugar_ring"_"$ru".txt
# Assign variable to this file
convergence=`cat convergence_"$type_sugar"_"$sugar_ring"_"$ru".txt`
fi

# Create variable that contains message that indicates that the calculation finished properly
string="    -- Stationary point found."

# Print NImag frequencies number
if [ "$sugar_ring" == "none" ]; then

grep -E -o 'NImag.{2}' "$type_sugar"_"$ru"_"$IL"_min.log > imag_freq_n.txt
# Count number of lines in file
cat imag_freq_n.txt | wc -l > number_lines_imag_freq_file.txt
lines_imag_file=`cat number_lines_imag_freq_file.txt`
if [ "$lines_imag_file" == 2 ]; then
sed -i 1,1d imag_freq_n.txt
fi

else

grep -E -o 'NImag.{2}' "$type_sugar"_"$sugar_ring"_"$ru"_"$IL"_min.log > imag_freq_n.txt
# Count number of lines in file
cat imag_freq_n.txt | wc -l > number_lines_imag_freq_file.txt
lines_imag_file=`cat number_lines_imag_freq_file.txt`
if [ "$lines_imag_file" == 2 ]; then
sed -i 1,1d imag_freq_n.txt
fi

fi

# Assign variable to the frequency file
freq_number=`cat imag_freq_n.txt`

# Conditions to determine convergence
if [ "$convergence" == "    -- Stationary point found." ];

then

if [ "$sugar_ring" == "none" ]; then
echo ""$environment"/"$type_sugar"/"$IL"/"$ru" converges" >> $MD_RUN/nucleotides/final_gaussian_opt/calc_that_did_converged_"$environment".log
else
echo ""$environment"/"$type_sugar"/"$sugar_ring"/"$IL"/"$ru" converges" >> $MD_RUN/nucleosides/final_gaussian_opt/calc_that_did_converged_"$environment".log
fi

# Determine if the optimized geometry is a min or TS of the PES
if [ "$freq_number" == "NImag=0" ]; then

echo "TC "$type_sugar" with "$sugar_ring" sugar ring conformation, "$IL" and RU "$ru" in "$environment" has NO Imaginary freq ($freq_number)" >> $MD_RUN/nucleosides/final_gaussian_opt/calc_with_no_imag_freq.log  

else # The case that there are imaginary frequencies

if [ "$sugar_ring" == "none" ]; then 
echo ""$environment"/"$type_sugar"/"$IL"/"$ru" has imaginary freq, #imag freq=($freq_number)" >> $MD_RUN/nucleosides/final_gaussian_opt/calc_with_imag_freq_"$environment".log  
else
echo ""$environment"/"$type_sugar"/"$sugar_ring"/"$IL"/"$ru" has imaginary freq, #imag freq=($freq_number)" >> $MD_RUN/nucleosides/final_gaussian_opt/calc_with_imag_freq_"$environment".log
fi

fi





else # THIS IS THE CASE THAT THE CALCULATION DID NOT CONVERGED






if [ "$sugar_ring" == "none" ]; then
echo ""$type_sugar"_"$ru".log in "$environment" did not converges" >> $MD_RUN/nucleosides/final_gaussian_opt/calc_that_did_not_converged_"$environment".log
# If the calculation does not converges then it is necesary to Convert log file to gjf
obabel "$type_sugar"_"$ru"_min.log -O nucleotide_nuevo.gjf
else
echo ""$type_sugar"_"$sugar_ring"_"$ru".log in "$environment" did not converges" >> $MD_RUN/nucleosides/final_gaussian_opt/calc_that_did_not_converged_"$environment".log
# If the calculation does not converges then it is necesary to Convert log file to gjf
obabel "$type_sugar"_"$sugar_ring"_"$ru"_min.log -O nucleotide_nuevo.gjf
fi

# Delete first two lines of the gjf
sed -i '1,2d' nucleoside_nuevo.gjf

# Paste instructions for b3lyp
cat $MD_RUN/instructions_DFT/instructions_b3lyp.txt nucleoside_nuevo.gjf > nucleoside_2opt.gjf

# Modify instructions if the modeling is in water
if [ "$environment" == "water" ]; then
sed -i 's/scf=maxcycle=1600/scrf=(solvent=water) scf=maxcycle=1600/g' nucleoside_2opt.gjf
fi

# Replace name chk by nucleoside name and new nucleoside name in job1.sh
if [ "$sugar_ring" == "none" ]; then
# Replace chk file in new gjf 2min file
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$type_sugar'_'$ru'_min.chk/g' nucleoside_2opt.gjf
# Replace name in job1.sh file
sed -i 's/g16 '$type_sugar'_'$ru'_min.gjf/g16 '$type_sugar'_'$ru'_2min.gjf/g' job1.sh
else
# Replace chk file in new gjf 2min file
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$type_sugar'_'$sugar_ring'_'$ru'_min.chk/g' nucleoside_2opt.gjf
# Replace name in job1.sh file
sed -i 's/g16 '$type_sugar'_'$sugar_ring'_'$ru'_min.gjf/g16 '$type_sugar'_'$sugar_ring'_'$ru'_2min.gjf/g' job1.sh
fi



# Add empty line at the end
echo "" >> nucleoside_2opt.gjf

# Change name of nucleoside
if [ "$sugar_ring" == "none" ]; then
mv nucleoside_2opt.gjf "$type_sugar"_"$ru"_2min.gjf
else
mv nucleoside_2opt.gjf "$type_sugar"_"$sugar_ring"_"$ru"_2min.gjf
fi





fi

}






##########################################################################################################################
##########################################################################################################################

############################################## WORK ON EACH FOLDER #######################################################

##########################################################################################################################
############################################## WORK ON EACH FOLDER #######################################################
##########################################################################################################################

# Open corresponding folder
cd nucleotides_a+b/final_gaussian_opt

for environment in "${environment[@]}"; do

mkdir $environment

cd $environment

for TC in "${TC_names[@]}"; do

cd $TC

if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

for ring_conf_fura in "${Ring_Conf_furanose[@]}"; do

cd $ring_conf_fura

for IL in "${IL_names[@]}"; do

cd $IL

for RU in "${RU_names[@]}"; do

cd $RU

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, ring conf $ring_conf_fura, IL $IL and RU $RU"

# Run script to obtain calculation status
completion_status $environment $TC $ring_conf_fura $IL $RU








cd ../
done

cd ../
done

cd ../
done

elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then

for ring_conf_pyra in "${Ring_Conf_pyranose[@]}"; do

cd $ring_conf_pyra

for IL in "${IL_names[@]}"; do

cd $IL

for RU in "${RU_names[@]}"; do

cd $RU

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, ring conf $ring_conf_pyra, IL $IL and RU $RU"

# Run script to obtain calculation status
completion_status $environment $TC $ring_conf_pyra $IL $RU










cd ../
done

cd ../
done

cd ../
done

elif [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then

for IL in "${IL_names[@]}"; do

cd $IL

for RU in "${RU_names[@]}"; do

cd $RU

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, none ring conf, IL $IL and RU $RU"

ring_conf="none"

# Run script to obtain calculation status
completion_status $environment $TC $ring_conf $IL $RU













cd ../
done

cd ../
done

fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT
