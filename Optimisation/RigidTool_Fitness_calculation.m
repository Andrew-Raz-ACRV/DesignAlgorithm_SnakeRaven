% Fitness Function Test script
% Andrew Razjigaev 2020 rigid tool
close all
clear all
clc
%% Anatomy Voxelization:
Anatomyfilename = 'VoxelDataMultiTarget.mat';
%Anatomyfilename = 'VoxelDataTask2Target1.mat';
%Anatomyfilename = 'VoxelDataAllTargets.mat';

%% Dexterity analysis settings:
sample_size = 10000000; 

%% Rigid Tool
rigid_results = zeros(1,10);
for ii=1:10
    %Rigid Tool
    disp(['Starting Fitness function test ' num2str(ii)])
    disp('Testing Design: Rigid Tool')

    tic
    Score_rigid = FastFitnessFunctionRigidTool(sample_size,Anatomyfilename);
    toc

    disp('Dexterity score for This design was:')
    disp(Score_rigid)
    rigid_results(ii) = Score_rigid;
end
disp('Dexterity scores for all Rigid Tool Analyses:')
disp(rigid_results)
save('DOF_Rigid_tool_analysis/Rigid_tool_baseline','rigid_results');
