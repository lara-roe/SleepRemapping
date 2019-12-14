#! /bin/bash
#SBATCH --job-name="fmriprep"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00

singularity run --cleanenv ~/my_images/fmriprep-1.4.0.simg ~/rds/rds-tb419-bekinschtein/Lara/SleepRemapping/ ~/rds/rds-tb419-bekinschtein/Lara/SleepRemapping/derivatives participant -w ~/rds/rds-tb419-bekinschtein/Lara/SleepRemapping/work/  -v --participant-label sub-S29 --fs-license-file ~/FreeSurfer/license.txt --output-space T1w 

