#!/bin/bash

#SBATCH --qos=priority
#SBATCH --job-name=Text
#SBATCH --account=impactee
#SBATCH --output=temp/temp-%j.out
#SBATCH --error=temp/temp-%j.err
#SBATCH --mail-user=maxkotz@pik-potsdam.de
#SBATCH --mail-type=fail
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00

#srun -n 1 python -u calc_low_clim_extremes.py tas
srun -n 1 python -u attribute_climate.py pr_extrm_95_num

