#!/bin/bash

#SBATCH --mem=16G             # memory, roughly 2 times %mem defined in the input name.com file
#SBATCH --time=01-00:00       # expect run time (DD-HH:MM)
#SBATCH --cpus-per-task=16    # No. of cpus for the job as defined by %nprocs in the name.com file
cat list_nucleosides.txt | while read nucleosides; 
do
module load gaussian/g16.b01
G16 $nucleosides           # G16 command, input: name.com, output: name.log
done
