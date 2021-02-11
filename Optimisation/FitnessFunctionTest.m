% Fitness Function Test script
% Andrew Razjigaev 2020 rigid tool
close all
clear all
clc
%% Anatomy Voxelization:
Anatomyfilename = 'VoxelDataMultiTarget.mat';
%Anatomyfilename = 'VoxelDataTask2Target1.mat';
%Anatomyfilename = 'VoxelDataAllTargets.mat';

%% Design Parameters:
% design.alpha = 0.65; %1.24;
% design.d = 3.49; %1.62;
% design.n = 10; %3;
% design.w = 4;
%0.3, 2,6,4 : 0.5, 1, 7, 4

%% Dexterity analysis settings:
sample_size = 10000000; 
% plotting = true;

% %% Fitness Function for Snake robot
% 
% disp('Starting Fitness function test . . .')
% disp(['Testing a ',num2str(length(design.alpha)),' Segment Design:'])
% disp(design)
% 
% %Function call debugging with figures:
% Score = FitnessFunctionVariableSegmentSnakeRobot(design,sample_size,plotting,Anatomyfilename);

% %Tell me number of available clusters:
% disp('Available Clusters:')
% ALLPROFILES = parallel.clusterProfiles;
% disp(ALLPROFILES)
% 
% %Number of Workers in the parallel pool:
% preferredNumWorkers = 30;
% 
% %Create the parallel pool of workers to do the for loop
% poolobj = parpool(preferredNumWorkers);
% 
% disp(poolobj)
% disp('Number of workers:')
% disp(preferredNumWorkers)
% 
% %Fitness function with parallel loop:
% tic
% Score = FastFitnessFunctionVariableSegmentSnakeRobot(design,sample_size,Anatomyfilename);
% toc
% 
% %End parallel loop delete the pool object
% delete(poolobj)
% 
% disp('Fitness score for This design was:')
% disp(Score)

% %% Fitness Function Results:
% result_file = strcat('Design_alpha',strrep(strrep(num2str(design.alpha),'.',''),' ','_'),...
%     '_n',strrep(strrep(num2str(design.n),'.',''),' ','_'),...
%     '_d',strrep(strrep(num2str(design.d),'.',''),' ','_'));
% 
% PlotResultsForFitnessFunction(result_file,Anatomyfilename)
% 
%% Rigid Tool
%Rigid Tool
disp('Starting Fitness function test . . .')
disp('Testing Design: Rigid Tool')
%Score_rigid = FitnessFunctionRigidTool(sample_size,plotting,Anatomyfilename);

tic
Score_rigid = FastFitnessFunctionRigidTool(sample_size,Anatomyfilename);
toc

disp('Dexterity score for This design was:')
disp(Score_rigid)

%% Fitness Function Results:
%result_file_rigid = 'Design_RigidTool_analysis_MultiTarget';
result_file_rigid = 'Design_RigidTool_analysis';

PlotResultsForFitnessFunction(result_file_rigid,Anatomyfilename)
