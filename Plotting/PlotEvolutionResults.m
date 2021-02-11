% Plot the results for Snake Evolution:
% Andrew Razjigaev February 2020

clc;
clear;
close all;

%% Load Data here:
Evolution_data = 'EvolutionResults07Apr2020.mat'; %Edit here
Evolution_data2 = 'EvolutionResults08Apr2020.mat'; %Edit here
V = load(Evolution_data);
V2 = load(Evolution_data2);

%Anatomy
Anatomyfilename = 'VoxelDataMultiTarget.mat';
% Anatomyfilename = 'VoxelDataTarget1_top.mat';
% Anatomyfilename = 'VoxelDataTarget2_right.mat';
% Anatomyfilename = 'VoxelDataTarget3_bottom.mat';
% Anatomyfilename = 'VoxelDataTarget4_left.mat';

%% Evolution Plot Results:

disp('Evolution Results. The Best solution was:');
disp(vector2designstruct(V.BestSol.Position))
disp('With best Dexterity:')
disp(-1*V.BestSol.Cost)
disp('\n');

%Plot Line of Evolution Best cost over time
figure;
plot(-1*V.BestCost,'LineWidth', 2); %semilogy(-1*V.BestCost, 'LineWidth', 2);
hold on
plot(-1*V2.BestCost,'LineWidth', 2); %semilogy(-1*V.BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

%% Best Solution Fitness Function Results

%Now get the best solution:
design = vector2designstruct(V.BestSol.Position);

%Get result file
result_file = strcat('Design_alpha',...
    strrep(strrep(strrep(strrep(num2str(design.alpha),'.',''),'        ','_'),'___','_'),'__','_'),...
    '_n',strrep(strrep(num2str(design.n),'.',''),'  ','_'),...
    '_d',strrep(strrep(strrep(strrep(num2str(design.d),'.',''),'        ','_'),'___','_'),'__','_'));


%Plot Visualisation of design, maximum service sphere and dexterity
%distribution
PlotResultsForFitnessFunction(result_file,Anatomyfilename)

