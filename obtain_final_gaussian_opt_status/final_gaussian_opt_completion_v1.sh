#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

function completion_status {

########################################################## Determine if the calculation converged or not #######################################################

# Print number of line of last instance of Maximum Force
if [ "$sugar_ring" == "none" ]; then
awk '/Maximum Force/ { ln = FNR } END { print ln }' "$type_sugar"_"$ru"_min.log > line_maximumforce_"$type_sugar"_"$ru".txt
# Assign a variable to this file 
line=`cat line_maximumforce_"$type_sugar"_"$ru".txt`
else
awk '/Maximum Force/ { ln = FNR } END { print ln }' "$type_sugar"_"$sugar_ring"_"$ru"_min.log > line_maximumforce_"$type_sugar"_"$sugar_ring"_"$ru".txt
# Assign a variable to this file 
line=`cat line_maximumforce_"$type_sugar"_"$sugar_ring"_"$ru".txt`
fi

# Obtain line that contains convergence message
total_lines=$(($line+6)) 

# Copy line with convergence message in new file
if [ "$sugar_ring" == "none" ]; then
sed -n ''$total_lines'p' "$type_sugar"_"$ru"_min.log > convergence_"$type_sugar"_"$ru".txt
# Assign variable to this file
convergence=`cat convergence_"$type_sugar"_"$ru".txt`
else
sed -n ''$total_lines'p' "$type_sugar"_"$sugar_ring"_"$ru"_min.log > convergence_"$type_sugar"_"$sugar_ring"_"$ru".txt
# Assign variable to this file
convergence=`cat convergence_"$type_sugar"_"$sugar_ring"_"$ru".txt`
fi

# Create variable that contains message that indicates that the calculation finished properly
string="    -- Stationary point found."

# Print NImag frequencies number
if [ "$sugar_ring" == "none" ]; then

grep -E -o 'NImag.{2}' "$type_sugar"_"$ru"_min.log > imag_freq_n.txt
# Count number of lines in file
cat imag_freq_n.txt | wc -l > number_lines_imag_freq_file.txt
lines_imag_file=`cat number_lines_imag_freq_file.txt`
if [ "$lines_imag_file" == 2 ]; then
sed -i 1,1d imag_freq_n.txt
fi

else

grep -E -o 'NImag.{2}' "$type_sugar"_"$sugar_ring"_"$ru"_min.log > imag_freq_n.txt
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
echo ""$type_sugar"_"$ru"_min.log in "$environment" converges" >> $MD_RUN/nucleosides/final_gaussian_opt/calc_that_did_converged.log
else
echo ""$type_sugar"_"$sugar_ring"_"$ru"_min.log in "$environment" converges" >> $MD_RUN/nucleosides/final_gaussian_opt/calc_that_did_converged.log
fi

# Determine if the optimized geometry is a min or TS of the PES
if [ "$freq_number" == "NImag=0" ]; then

echo "opt + freq gaussian calc of TC "$type_sugar" with "$sugar_ring" sugar ring conformation and RU "$ru" in "$environment" has NO Imaginary freq ($freq_number)" >> $MD_RUN/nucleosides/final_gaussian_opt/calc_with_no_imag_freq.log  

else # The case that there are imaginary frequencies

echo "opt + freq gaussian calc of TC "$type_sugar" with "$sugar_ring" sugar ring conformation and RU "$ru" in "$environment" has imaginary freq ($freq_number)" >> $MD_RUN/nucleosides/final_gaussian_opt/calc_with_imag_freq.log  

fi





else # THIS IS THE CASE THAT THE CALCULATION DID NOT CONVERGED






if [ "$sugar_ring" == "none" ]; then
echo ""$type_sugar"_"$ru".log in "$environment" did not converges" >> $MD_RUN/nucleosides/final_gaussian_opt/calc_that_did_not_converged.log
# If the calculation does not converges then it is necesary to Convert log file to gjf
obabel "$type_sugar"_"$ru"_min.log -O nucleoside_nuevo.gjf
else
echo ""$type_sugar"_"$sugar_ring"_"$ru".log in "$environment" did not converges" >> $MD_RUN/nucleosides/final_gaussian_opt/calc_that_did_not_converged.log
# If the calculation does not converges then it is necesary to Convert log file to gjf
obabel "$type_sugar"_"$sugar_ring"_"$ru"_min.log -O nucleoside_nuevo.gjf
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

# Open corresponding folder
cd nucleosides/final_gaussian_opt

echo "vacuum" > environment.txt
echo "water" >> environment.txt

cat environment.txt | while read environment

do

cd $environment

cat $MD_RUN/list_folders/sugars_list.txt | while read sugars

do

if [[ "$sugars" == *\/* ]]; then # This is the case that the TC is a sugar with sugar ring conformation

# Obtain trifunctional connector name
echo $sugars | cut -d "/" -f1 > sugar_name.txt
type_sugar=`cat sugar_name.txt`

# Obtain sugar ring conformation name
echo $sugars | cut -d "/" -f2 > sugar_ring_conf.txt
sugar_ring=`cat sugar_ring_conf.txt`

cd $type_sugar/$sugar_ring

cat $MD_RUN/list_folders/bases_list.txt | while read ru

do

cd $ru

# Message
echo "working on $environment environment, $type_sugar with $sugar_ring sugar ring conformation and base $ru"

# Run function that determines if calculation completed
completion_status






cd ../

done

cd ../../

else # This is the case that the TC is glycerol or glyceric acid or peptide which is not a sugar and does not has sugar ring conformation

# Obtain trifunctional connector name
echo $sugars > sugar_name.txt
type_sugar=`cat sugar_name.txt`

# Obtain sugar ring conformation name
sugar_ring="none"

if [ "$type_sugar" != "peptide" ]; then

cd $type_sugar

cat $MD_RUN/list_folders/bases_list.txt | while read ru

do

cd $ru

sugar_ring="none"

# Message
echo "working on $environment environment, $type_sugar with $sugar_ring sugar ring conformation and base $ru"

# Run function that determines if calculation completed
completion_status






cd ../

done

elif [ "$type_sugar" == "peptide" ]; then

cd $type_sugar

cat $MD_RUN/list_folders/bases_ac_list.txt | while read ru

do

cd $ru

sugar_ring="none"

# Message
echo "working on $environment environment, $type_sugar with $sugar_ring sugar ring conformation and base $ru"

# Run function that determines if calculation completed
completion_status









cd ../

done

fi


cd ../

fi

done

cd ../

done
