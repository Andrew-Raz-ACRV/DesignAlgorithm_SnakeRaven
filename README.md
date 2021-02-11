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
You can modify the settings for different task objectives

## 4) Results:
