# DesignAlgorithm_SnakeRaven
Automated design algorithms for end-to-end design optimisation of a patient-specific surgical manipulator considering dexterity, reachability and task-space obstacles. 

## Prerequisites:
You will need MATLAB and its parallel computing toolbox

## 1) Patient Scan
You will need an stl file. The one we used can be downloaded through this link:

## 2) Voxelisation
You will need to run this file to generate the voxel map:
```
GenerateVoxelisation.m
```
The output of this file 

## 3) Optimisation:
The main file for optimisation is:
```
SnakeRaven_Evolution_script.m
```
You can run this on a HPC or another powerful computer utilising a parpool.
You can modify the settings for different task objectives by commenting out different options
The file will create a new directory which saves all fitness function result files e.g. Design_alphaXXX_XXX_nXX_XX_dXXX_XXX.mat
Where XXX fills ths pace of the actual design parameters
And the final overall results: Snake_Evolution_ResultsXX-XXX_XXXX_XX_XX_XX.mat
Where XXX fills the date and time of the complete evolution

## 4) Results:
