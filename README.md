# DesignAlgorithm_SnakeRaven
Automated design algorithms for end-to-end design optimisation of a patient-specific surgical manipulator considering dexterity, reachability and task-space obstacles. Based on our publication "End-to-end design of a bespoke 3D printed dexterous snake-like surgical manipulator using an evolution algorithm"

## Prerequisites:
You will need MATLAB and its parallel computing toolbox

## 1) Patient Scan
You will need an stl file of the anatomy. The one we used can be downloaded through this link: https://www.dropbox.com/s/bmo72qa6k2d3a8v/kneemodel.stl?dl=0

## 2) Voxelisation
You will need to run this file to generate the voxel map:
```
GenerateVoxelisation.m
```
An example of the output of voxelisation is 'VoxelDataMultiTarget.mat' but it will have the date in the name instead

## 3) Optimisation:
The main file for optimisation is:
```
SnakeRaven_Evolution_script.m
```
You can run this on a high-performance cluster (HPC) or another powerful computer utilising a parpool.
You can modify the settings for different task objectives by commenting out different options e.g. 

```
%Anatomyfilename = 'VoxelDataMultiTarget.mat'; 
```

The code will create a new directory that saves all fitness function result files e.g. Design_alphaXXX_XXX_nXX_XX_dXXX_XXX.mat
Where XXX fills the space of the actual design parameters

The final overall results are saved as a file called: Snake_Evolution_ResultsXX-XXX_XXXX_XX_XX_XX.mat
Where XXX fills in the date and time of the complete evolution. 
It will also create a backup results file after each generation.

The backup results can be used to continue executing the evolution algorithm if it fails remotely.
You will need to edit it (so that it knows where the directory is) and then run the script:
```
Revive_Evolution.m
```

## Example of using this script on HPC in a PBS job system:

To run the code on HPC you can create a .sh file and submit the job to the cluster. 
Below is an example of a .sh file let's call it "a_pbs_job.sh" 
This script will open Matlab 2018b and run the evolution algorithm on HPC using 1 core, 12 parallel processors, 5GB of memory and an execution time limit of 200 hours.
```
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
```

To execute it, log in to HPC and navigate to your directory.
Edit the PBS script and for thoroughness ensure the text is compatible with Linux by running: 
```
dos2unix a_pbs_job.sh
```
Afterwards, you can submit the job request via the command line:
```
$ qsub a_pbs_job.sh
```
Which returns an ID for that particular job.
You can check the status of all your jobs via your username:
```
$ qstat -USERNAME
```
You can abort the particular job using the job ID:
```
$ qdel job_id
```

## 4) Plotting:
You can now plot all the evolution results to get:
- a design render
- fitness over time
- mean fitness and standard deviation over time
- boxplot of the parameter variation
- dexterity distribution graph
- maximum service sphere

see the example code for how to plot the results 
```
PlotEvolutionResultsSnake_Example.m
```

## Patient Specific Design Process:
![alt text](https://github.com/Andrew-Raz-ACRV/DesignAlgorithm_SnakeRaven/blob/main/Plotting/Patient_specific_Flowchart_pictures-2.jpg)
