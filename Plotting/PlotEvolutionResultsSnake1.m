% Plot the results for the DOF One Segment Snake Evolution Experiments:
% Andrew Razjigaev 24 April 2020

clc;
clear;
close all;

%% Load Data here:
Snake1 = cell(10,1);
Snake1{1} = 'Snake_Evolution_Results16-Apr-2020_12_53_04 Finished_17-Apr-2020_23_55_39.mat';
Snake1{2} = 'Snake_Evolution_Results16-Apr-2020_12_53_14 Finished_18-Apr-2020_02_38_28.mat';
Snake1{3} = 'Snake_Evolution_Results16-Apr-2020_12_53_23 Finished_18-Apr-2020_19_09_32.mat';
Snake1{4} = 'Snake_Evolution_Results16-Apr-2020_13_14_37 Finished_18-Apr-2020_20_10_03.mat';
Snake1{5} = 'Snake_Evolution_Results16-Apr-2020_13_14_40 Finished_18-Apr-2020_02_26_20.mat';
Snake1{6} = 'Snake_Evolution_Results16-Apr-2020_13_35_33 Finished_17-Apr-2020_19_54_13.mat';
Snake1{7} = 'Snake_Evolution_Results16-Apr-2020_13_44_58 Finished_18-Apr-2020_03_46_44.mat';
Snake1{8} = 'Snake_Evolution_Results16-Apr-2020_13_44_54 Finished_18-Apr-2020_02_49_46.mat';
Snake1{9} = 'Snake_Evolution_Results16-Apr-2020_13_45_09 Finished_18-Apr-2020_02_59_26.mat';
Snake1{10} = 'Snake_Evolution_Results16-Apr-2020_14_15_21 Finished_18-Apr-2020_03_42_50.mat';

% Primary Dataset
V = load(Snake1{1});
directory = 'Snake_Evolution_Results16-Apr-2020_12_53_04';

%Anatomy
Anatomyfilename = 'VoxelDataMultiTarget.mat';
% Anatomyfilename = 'VoxelDataTarget1_top.mat';
% Anatomyfilename = 'VoxelDataTarget2_right.mat';
% Anatomyfilename = 'VoxelDataTarget3_bottom.mat';
% Anatomyfilename = 'VoxelDataTarget4_left.mat';

%% Evolution Fitness over time Plot:

disp('Evolution Results. The Best solution was:');
disp(vector2designstruct(V.BestSol.Position))
disp('With best Dexterity:')
disp(-1*V.BestSol.Cost)
disp('\n');

%Plot Line of Evolution Best cost over time
figure;
for ii = 1:10
    Data = load(Snake1{ii});
    plot(-1*Data.BestCost,'LineWidth', 2);
    hold on
end
title('Fitness Over Time');
xlabel('Iteration');
ylabel('Best Cost');
grid on;

%Plot Mean and standard deviation range Lines of Evolution Best cost over time
All_data = zeros(100,10);
for ii = 1:10
    Data = load(Snake1{ii});
    All_data(:,ii) = -1*Data.BestCost;
end

mean_fitness = mean(All_data,2);
sigma_fitness = std(All_data,0,2);
Upper_line = mean_fitness + sigma_fitness;
Lower_line = mean_fitness - sigma_fitness;

figure;
plot(Upper_line,'b-','LineWidth', 2); hold on;
plot(Lower_line,'b-','LineWidth', 2); hold on;
%plot shaded area
plot(mean_fitness,'r-','LineWidth', 2); hold on;
title('Mean Fitness and Standard Deviation Over Time');
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
PlotResultsForFitnessFunction(strcat(directory,'/',result_file),Anatomyfilename)

