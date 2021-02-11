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
The output of this file like 'VoxelDataMultiTarget.mat' but it will have the date in the name instead

## 3) Optimisation:
The main file for optimisation is:
```
SnakeRaven_Evolution_script.m
```
You can run this on a HPC or another powerful computer utilising a parpool.
You can modify the settings for different task objectives by commenting out different options e.g. 

```
%Anatomyfilename = 'VoxelDataMultiTarget.mat'; 
```

The code will create a new directory which saves all fitness function result files e.g. Design_alphaXXX_XXX_nXX_XX_dXXX_XXX.mat
Where XXX fills the space of the actual design parameters

The final overall results: Snake_Evolution_ResultsXX-XXX_XXXX_XX_XX_XX.mat
Where XXX fills the date and time of the complete evolution. 
It will also create a back up results file after each generation.

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
