#!/bin/bash -l

#### Job Name
#PBS -N SnakeRaven

####Select Resources

#PBS -l nodes=1:ppn=12
#PBS -l mem=5gb
#PBS -l walltime=200:00:00

#### Load Matlab module (setup enviornment)
module load matlab/2018b

#### Execute program
matlab -r SnakeRaven_Evolution_script -logfile logfile_SnakeRaven_Evolution_script.log -nodisplay -nodesktop -nosplash