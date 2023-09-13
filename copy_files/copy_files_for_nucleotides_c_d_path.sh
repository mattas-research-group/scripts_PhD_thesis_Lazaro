#!/bin/bash

# Obtain current location of this script
MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# List of environments, TC, ring conformations and RUs
declare -a environment=("vacuum" "water")
declare -a TC_names=("2deoxy_ribofuranose" "ribofuranose" "2deoxy_ribopyranose" "ribopyranose" "threose" "glycerol" "glyceric_acid")
declare -a Ring_Conf_furanose=("2endo_alpha" "2endo_beta" "3endo_alpha" "3endo_beta")
declare -a Ring_Conf_pyranose=("D1C4_alpha" "D1C4_beta" "D4C1_alpha" "D4C1_beta")
declare -a RU_names=("adenine" "guanine" "cytosine" "thymine" "uracil" "TARC_cbond" "TARC_nbond" "BA_cbond" "BA_nbond" "CA" "melamine")
declare -a IL_names=("phosphate" "arsenate")


# Generate function that copies corresponding files
function copy_files {

# Load corresponding variables
local environment=$1
local TC=$2
local ring_conf=$3
local IL=$4
local RU=$5

# Copy RUs final gaussian opt log file
if [ $RU == "BA_cbond" ] || [ "$RU" == "BA_nbond" ] || [ "$RU" == "TARC_cbond" ] || [ "$RU" == "TARC_nbond" ] || [ "$RU" == "melamine" ]; then
cp $MD_RUN/folders_to_copy/RUs/final_gaussian_opt/$environment/$RU/$RU.mol . 
else
cp $MD_RUN/folders_to_copy/RUs/final_gaussian_opt/$environment/$RU/$RU.log .
fi

# Load nucleoside final gaussian opt log file

if [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then

# Load the files
backbone_1opt=$MD_RUN/TCs_ILs_backbone/final_gaussian_opt/$environment/$TC/$IL/"$TC"_"$IL"_min.log
backbone_2opt=$MD_RUN/TCs_ILs_backbone/final_gaussian_opt/$environment/$TC/$IL/"$TC"_"$IL"_min2.log
backbone_3opt=$MD_RUN/TCs_ILs_backbone/final_gaussian_opt/$environment/$TC/$IL/"$TC"_"$IL"_min3.log

else

# Load the files
backbone_1opt=$MD_RUN/TCs_ILs_backbone/final_gaussian_opt/$environment/$TC/$ring_conf/$IL/"$TC"_"$ring_conf"_"$IL"_min.log
backbone_2opt=$MD_RUN/TCs_ILs_backbone/final_gaussian_opt/$environment/$TC/$ring_conf/$IL/"$TC"_"$ring_conf"_"$IL"_min2.log
backbone_3opt=$MD_RUN/TCs_ILs_backbone/final_gaussian_opt/$environment/$TC/$ring_conf/$IL/"$TC"_"$ring_conf"_"$IL"_min3.log

fi


# Copy RUs-ILs backbone final gaussian opt log file:

if [ -f "$backbone_3opt" ]; then

log_file=$backbone_3opt

echo "there was a 3nd opt"

if [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then
echo ""$environment"/"$TC"_"$IL" ------------------> min3" >> $MD_RUN/number_runs.txt
else
echo ""$environment"/"$TC"_"$ring_conf"_"$IL" ------------------> min3" >> $MD_RUN/number_runs.txt
fi

else

if [ -f "$backbone_2opt" ]; then

log_file=$backbone_2opt

echo "there was a 2nd opt"

if [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then
echo ""$environment"/"$TC"_"$IL" ------------------> min2" >> $MD_RUN/number_runs.txt
else
echo ""$environment"/"$TC"_"$ring_conf"_"$IL" ------------------> min2" >> $MD_RUN/number_runs.txt
fi

else

log_file=$backbone_1opt

echo "there was only 1 opt" 

if [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then
echo ""$environment"/"$TC"_"$IL" ------------------> min" >> $MD_RUN/number_runs.txt
else
echo ""$environment"/"$TC"_"$ring_conf"_"$IL" ------------------> min" >> $MD_RUN/number_runs.txt
fi


fi
fi

}




############################################## WORK ON EACH FOLDER #######################################################

# Open corresponding folder
mkdir nucleotides_c+d
mkdir nucleotides_c+d/gaussian_scan
cd nucleotides_c+d/gaussian_scan

for environment in "${environment[@]}"; do

mkdir $environment

cd $environment

for TC in "${TC_names[@]}"; do

mkdir $TC
cd $TC

if [ "$TC" == "2deoxy_ribofuranose" ] || [ "$TC" == "ribofuranose" ] || [ "$TC" == "threose" ]; then

for ring_conf_fura in "${Ring_Conf_furanose[@]}"; do

mkdir $ring_conf_fura

cd $ring_conf_fura

for IL in "${IL_names[@]}"; do

mkdir $IL
cp $log_file .
cp $log_file "$TC"_"$ring_conf_fura"_"$RU".log
cd $IL

for RU in "${RU_names[@]}"; do

mkdir $RU

cd $RU

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, ring conf $ring_conf_fura, IL $IL and RU $RU"

# Time to copy the corresponding files from the nucleosides folders and the IL folders
copy_files $environment $TC $ring_conf_fura $IL $RU 

cp $log_file .
cp $log_file "$TC"_"$ring_conf_fura"_"$IL".log




cd ../
done

cd ../
done

cd ../
done

elif [ "$TC" == "2deoxy_ribopyranose" ] || [ "$TC" == "ribopyranose" ]; then

for ring_conf_pyra in "${Ring_Conf_pyranose[@]}"; do

mkdir $ring_conf_pyra

cd $ring_conf_pyra

for IL in "${IL_names[@]}"; do

mkdir $IL

cd $IL

for RU in "${RU_names[@]}"; do

mkdir $RU

cd $RU

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, ring conf $ring_conf_pyra, IL $IL and RU $RU"

# Time to copy the corresponding files from the nucleosides folders and the IL folders
copy_files $environment $TC $ring_conf_pyra $IL $RU 

cp $log_file .
cp $log_file "$TC"_"$ring_conf_pyra"_"$IL".log





cd ../
done

cd ../
done

cd ../
done

elif [ "$TC" == "glycerol" ] || [ "$TC" == "glyceric_acid" ]; then

for IL in "${IL_names[@]}"; do

mkdir $IL

cd $IL

for RU in "${RU_names[@]}"; do

mkdir $RU

cd $RU

# Print message of on where folder is been worked
echo "working on $environment, TC $TC, none ring conf, IL $IL and RU $RU"

ring_conf="none"

# Time to copy the corresponding files from the nucleosides folders and the IL folders
copy_files $environment $TC $ring_conf $IL $RU

cp $log_file .
cp $log_file "$TC"_"$IL".log



cd ../
done

cd ../
done





fi

cd ../ # LEAVE CURRENT TC FOLDER
done # END OF CYCLE THROUGH TC FOLDERS

cd ../ # LEAVE CURRENT ENVIRONMENT FOLDER
done  # END OF CYCLE THROUGH ENVIRONMENT
