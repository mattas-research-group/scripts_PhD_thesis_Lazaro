#!/bin/bash

#SBATCH --mem=16G             # memory, roughly 2 times %mem defined in the input name.com file
#SBATCH --time=00-05:00       # expect run time (DD-HH:MM)
#SBATCH --cpus-per-task=16    # No. of cpus for the job as defined by %nprocs in the name.com file
cd b3lyp
cat list_b3lyp.txt | while read n; 
do
module load gaussian/g16.b01
G16 $n.gjf            # G16 command, input: name.com, output: name.log
done
cd ../
cd m052x
cat list_m052x.txt | while read m; 
do
module load gaussian/g16.b01
G16 $m.gjf            # G16 command, input: name.com, output: name.log
done
cd ../
