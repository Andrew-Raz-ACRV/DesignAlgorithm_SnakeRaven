% Fitness Function Plot The robot
% Andrew Razjigaev 2020 Snake Robot
% This scripts runs the plotting version of the fitness function and
% exports the figure into high resolution!

close all
clear all
clc

%% Anatomy Voxelization:
Anatomyfilename = 'VoxelDataMultiTarget.mat'; %T1234
%Anatomyfilename = 'VoxelDataTarget1_top.mat'; %T1
%Anatomyfilename = 'VoxelDataTarget2_right.mat'; %T2
%Anatomyfilename = 'VoxelDataTarget3_bottom.mat'; %T3
%Anatomyfilename = 'VoxelDataTarget4_left.mat'; %T4
%Anatomyfilename = 'VoxelDataTask2Target1.mat'; %T5
Voxels = load(Anatomyfilename,'Voxel_data');

%% Designs:

% %One:
% design.alpha = 0.8700;
% design.n = 3;
% design.d = 1.1200;
% design.w = 4;
% 
%Two
design.alpha = [0.2000 0.8800];
design.n = [3 3];
design.d = [1 1];
design.w = [4 4];
% 
% %T1
% design.alpha = [0.1300 1.1300];
% design.n = [5 3];
% design.d = [1.3500 1.0600];
% design.w = [4 4];
% 
% %T2
% design.alpha = [0.0500 0.6700];
% design.n = [6 4];
% design.d = [2.8900 1];
% design.w = [4 4];
% 
% %T3
% design.alpha = [0.2600 1.0100];
% design.n = [1 2];
% design.d = [7.9400 1.0100];
% design.w = [4 4];
% 
% %T4
% design.alpha = [0.2900 1.1000];
% design.n = [1 3];
% design.d = [2.3800 1];
% design.w = [4 4];

%T5
% design.alpha = [0.9500 0.9900];
% design.n = [3 2];
% design.d = [9.1500 1.0400];
% design.w = [4 4];

%% Run Calculation

sample_size = 100000; 
directory = 'EVOLUTION_RESULTS_2020';

disp('Starting Fitness function test . . .')
disp(['Testing a ',num2str(length(design.alpha)),' Segment Design:'])
disp(design)

Score = PlotFitnessFunctionVariableSegmentSnakeRobot(design,sample_size,Voxels,directory);

%% Ends during the Fitness Function
% %% Fitness Function Results:
% result_file = strcat('Design_alpha',strrep(strrep(num2str(design.alpha),'.',''),' ','_'),...
%     '_n',strrep(strrep(num2str(design.n),'.',''),' ','_'),...
%     '_d',strrep(strrep(num2str(design.d),'.',''),' ','_'));
% 
% PlotResultsForFitnessFunction(result_file,Anatomyfilename)