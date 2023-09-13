#!/bin/bash
#SBATCH --mem=16G             # memory, roughly 2 times %mem defined in the input name.com file
#SBATCH --time=00-20:00       # expect run time (DD-HH:MM)
#SBATCH --cpus-per-task=16    # No. of cpus for the job as defined by %nprocs in the name.com file
module load gaussian/g16.b01
g16 2deoxy_ribofuranose_2endo_alpha_guanine_phosphate_scan2.gjf            # G16 command, input: name.com, output: name.log
wait
