#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# Create general folder for gaussian scans
mkdir $MD_RUN/nucleosides/gaussian_scan

############################################## WORK ON EACH FOLDER #######################################################

# Open corresponding folder
cd nucleosides/initial_gaussian_opt

echo "vacuum" > environment.txt
echo "water" >> environment.txt

cat environment.txt | while read environment

do

# Create environment folder in folder for scans
mkdir $MD_RUN/nucleosides/gaussian_scan/$environment

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

# Create sugar and sugar ring folders in folders for scan
mkdir $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar
mkdir $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/$sugar_ring


cd $type_sugar/$sugar_ring

cat $MD_RUN/list_folders/bases_list.txt | while read ru

do

mkdir $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/$sugar_ring/$ru


cd $ru

# Message
echo "working on $environment environment, $type_sugar with $sugar_ring sugar ring conformation and base $ru"

########################################################## Determine if the calculation converged or not #######################################################

# Print number of line of last instance of Maximum Force

awk '/Maximum Force/ { ln = FNR } END { print ln }' "$type_sugar"_"$sugar_ring"_"$ru".log > line_maximumforce_"$type_sugar"_"$sugar_ring"_"$ru".txt

# Assign a variable to this file 
line=`cat line_maximumforce_"$type_sugar"_"$sugar_ring"_"$ru".txt`

# Obtain line that contains convergence message
total_lines=$(($line+6)) 

# Copy line with convergence message in new file
sed -n ''$total_lines'p' "$type_sugar"_"$sugar_ring"_"$ru".log > convergence_"$type_sugar"_"$sugar_ring"_"$ru".txt
# Assign variable to this file
convergence=`cat convergence_"$type_sugar"_"$sugar_ring"_"$ru".txt`

# Create variable that contains message that indicates that the calculation finished properly
string="    -- Stationary point found."

# Conditions to determine convergence
if [ "$convergence" == "    -- Stationary point found." ];

then

echo ""$type_sugar"_"$sugar_ring"_"$ru".log in "$environment" converges" >> $MD_RUN/nucleosides/initial_gaussian_opt/convergence_results.log

# Convert this file to mol file and copy it in scan folders
obabel "$type_sugar"_"$sugar_ring"_"$ru".log -O $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/$sugar_ring/$ru/"$type_sugar"_"$sugar_ring"_"$ru".mol


else


echo ""$type_sugar"_"$sugar_ring"_"$ru".log in "$environment" DOES NOT CONVERGES" >> $MD_RUN/nucleosides/initial_gaussian_opt/convergence_results.log

# If the calculation does not converges then it is necesary to Convert log file to gjf
obabel "$type_sugar"_"$sugar_ring"_"$ru".log -O "$type_sugar"_"$sugar_ring"_"$ru"_nuevo.gjf

# Delete first two lines of the gjf
sed -i '1,3d' "$type_sugar"_"$sugar_ring"_"$ru"_nuevo.gjf

# Paste instructions for b3lyp
cat $MD_RUN/instructions_DFT/instructions_b3lyp.txt "$type_sugar"_"$sugar_ring"_"$ru"_nuevo.gjf > "$type_sugar"_"$sugar_ring"_"$ru"_2opt.gjf

# Modify instructions if the modeling is in water
if [ "$environment" == "water" ]; then
sed -i 's/scf=maxcycle=1600/scrf=(solvent=water) scf=maxcycle=1600/g' "$type_sugar"_"$sugar_ring"_"$ru"_2opt.gjf
fi

# Replace name chk by nucleoside name
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$type_sugar'_'$sugar_ring'_'$ru'.chk/g' "$type_sugar"_"$sugar_ring"_"$ru"_2opt.gjf

# Print this cell in new list
echo "$type_sugar"_"$sugar_ring"_"$ru"_2opt.gjf >> list_new_b3lyp.txt

# Add frozen coordinates
if [ "$type_sugar" == "2deoxy_ribofuranose" ] || [ "$type_sugar" == "ribofuranose" ]; then
cat "$type_sugar"_"$sugar_ring"_"$ru"_2opt.gjf $MD_RUN/frozen_coordinates/frozen_coord_fura.txt > "$type_sugar"_"$sugar_ring"_"$ru"_2opt_new.gjf
mv "$type_sugar"_"$sugar_ring"_"$ru"_2opt_new.gjf "$type_sugar"_"$sugar_ring"_"$ru"_2opt.gjf
elif [ "$type_sugar" == "threose" ]; then
cat "$type_sugar"_"$sugar_ring"_"$ru"_2opt.gjf $MD_RUN/frozen_coordinates/frozen_coord_threo.txt > "$type_sugar"_"$sugar_ring"_"$ru"_2opt_new.gjf
mv "$type_sugar"_"$sugar_ring"_"$ru"_2opt_new.gjf "$type_sugar"_"$sugar_ring"_"$ru"_2opt.gjf
elif [ "$type_sugar" == "2deoxy_ribopyranose" ] || [ "$type_sugar" == "ribopyranose" ]; then
cat "$type_sugar"_"$sugar_ring"_"$ru"_2opt.gjf $MD_RUN/frozen_coordinates/frozen_coord_pyra.txt > "$type_sugar"_"$sugar_ring"_"$ru"_2opt_new.gjf
mv "$type_sugar"_"$sugar_ring"_"$ru"_2opt_new.gjf "$type_sugar"_"$sugar_ring"_"$ru"_2opt.gjf
fi


fi






cd ../

done

cd ../../

else # This is the case that the TC is glycerol or glyceric acid or peptide which is not a sugar and does not has sugar ring conformation

# Obtain trifunctional connector name
echo $sugars > sugar_name.txt
type_sugar=`cat sugar_name.txt`

# Obtain sugar ring conformation name
sugar_ring="none"

# Create sugar and sugar ring folders in folders for scan
mkdir $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar


if [ "$type_sugar" != "peptide" ]; then

cd $type_sugar

cat $MD_RUN/list_folders/bases_list.txt | while read ru

do

# Create RU folder in gaussian scan folder
mkdir $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/$ru


cd $ru

sugar_ring="none"

# Message
echo "working on $environment environment, $type_sugar with $sugar_ring sugar ring conformation and base $ru"

########################################################## Determine if the calculation converged or not #######################################################

# Print number of line of last instance of Maximum Force

awk '/Maximum Force/ { ln = FNR } END { print ln }' "$type_sugar"_"$ru".log > line_maximumforce_"$type_sugar"_"$ru".txt

# Assign a variable to this file 
line=`cat line_maximumforce_"$type_sugar"_"$ru".txt`

# Obtain line that contains convergence message
total_lines=$(($line+6)) 

# Copy line with convergence message in new file
sed -n ''$total_lines'p' "$type_sugar"_"$ru".log > convergence_"$type_sugar"_"$ru".txt
# Assign variable to this file
convergence=`cat convergence_"$type_sugar"_"$ru".txt`

# Create variable that contains message that indicates that the calculation finished properly
string="    -- Stationary point found."

# Conditions to determine convergence
if [ "$convergence" == "    -- Stationary point found." ];

then

echo ""$type_sugar"_"$ru".log in "$environment" converges" >> $MD_RUN/nucleosides/initial_gaussian_opt/convergence_results.log

# Convert this file to mol file and copy it in scan folders
obabel "$type_sugar"_"$ru".log -O $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/$ru/"$type_sugar"_"$ru".mol


else


echo ""$type_sugar"_"$ru".log in "$environment" DOES NOT CONVERGES" >> $MD_RUN/nucleosides/initial_gaussian_opt/convergence_results.log

# If the calculation does not converges then it is necesary to Convert log file to gjf
obabel "$type_sugar"_"$ru".log -O "$type_sugar"_"$ru"_nuevo.gjf

# Delete first three lines of the gjf
sed -i '1,3d' "$type_sugar"_"$ru"_nuevo.gjf

# Paste instructions for b3lyp
cat $MD_RUN/instructions_DFT/instructions_b3lyp.txt "$type_sugar"_"$ru"_nuevo.gjf > "$type_sugar"_"$ru"_2opt.gjf

# Modify instructions if the modeling is in water
if [ "$environment" == "water" ]; then
sed -i 's/scf=maxcycle=1600/scrf=(solvent=water) scf=maxcycle=1600/g' "$type_sugar"_"$ru"_2opt.gjf
fi

# Replace name chk by nucleoside name
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$type_sugar'_'$ru'.chk/g' "$type_sugar"_"$ru"_2opt.gjf

# Replace keyword modredundant 
sed -i 's/opt=(readfc,modredundant)/opt=readfc/g' "$type_sugar"_"$ru"_2opt.gjf

# Print this cell in new list
echo "$type_sugar"_"$ru"_2opt.gjf >> list_new_b3lyp.txt


fi





cd ../

done

elif [ "$type_sugar" == "peptide" ]; then

cd $type_sugar

cat $MD_RUN/list_folders/bases_ac_list.txt | while read ru_ac

do

# Create RU folder in gaussian scan folder
mkdir $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/$ru_ac


cd $ru_ac

sugar_ring="none"

# Message
echo "working on $environment environment, $type_sugar with $sugar_ring sugar ring conformation and base $ru_ac"

########################################################## Determine if the calculation converged or not #######################################################

# Print number of line of last instance of Maximum Force

awk '/Maximum Force/ { ln = FNR } END { print ln }' "$type_sugar"_"$ru_ac".log > line_maximumforce_"$type_sugar"_"$ru_ac".txt

# Assign a variable to this file 
line=`cat line_maximumforce_"$type_sugar"_"$ru_ac".txt`

# Obtain line that contains convergence message
total_lines=$(($line+6)) 

# Copy line with convergence message in new file
sed -n ''$total_lines'p' "$type_sugar"_"$ru_ac".log > convergence_"$type_sugar"_"$ru_ac".txt
# Assign variable to this file
convergence=`cat convergence_"$type_sugar"_"$ru_ac".txt`

# Create variable that contains message that indicates that the calculation finished properly
string="    -- Stationary point found."

# Conditions to determine convergence
if [ "$convergence" == "    -- Stationary point found." ];

then

echo ""$type_sugar"_"$ru_ac".log in "$environment" converges" >> $MD_RUN/nucleosides/initial_gaussian_opt/convergence_results.log

# Convert this file to mol file and copy it in scan folders
obabel "$type_sugar"_"$ru_ac".log -O $MD_RUN/nucleosides/gaussian_scan/$environment/$type_sugar/$ru/"$type_sugar"_"$ru_ac".mol


else


echo ""$type_sugar"_"$ru_ac".log in "$environment" DOES NOT CONVERGES" >> $MD_RUN/nucleosides/initial_gaussian_opt/convergence_results.log

# If the calculation does not converges then it is necesary to Convert log file to gjf
obabel "$type_sugar"_"$ru_ac".log -O "$type_sugar"_"$ru_ac"_nuevo.gjf

# Delete first three lines of the gjf
sed -i '1,3d' "$type_sugar"_"$ru_ac"_nuevo.gjf

# Paste instructions for b3lyp
cat $MD_RUN/instructions_DFT/instructions_b3lyp.txt "$type_sugar"_"$ru_ac"_nuevo.gjf > "$type_sugar"_"$ru_ac"_2opt.gjf

# Modify instructions if the modeling is in water
if [ "$environment" == "water" ]; then
sed -i 's/scf=maxcycle=1600/scrf=(solvent=water) scf=maxcycle=1600/g' "$type_sugar"_"$ru_ac"_2opt.gjf
fi

# Replace name chk by nucleoside name
sed -i 's/%chk=2deoxyribo_beta_vacuum_cell21.chk/%chk='$type_sugar'_'$ru_ac'.chk/g' "$type_sugar"_"$ru_ac"_2opt.gjf

# Replace keyword modredundant 
sed -i 's/opt=(readfc,modredundant)/opt=readfc/g' "$type_sugar"_"$ru_ac"_2opt.gjf

# Print this cell in new list
echo "$type_sugar"_"$ru_ac"_2opt.gjf >> list_new_b3lyp.txt


fi






cd ../

done

fi


cd ../

fi

done

cd ../

done




