% Fitness Function tests for all generated designs for T5 or any other target:
% Andrew Razjigaev 2020 using all best designs

%Test Runs:
%Rigid Tool on T5 x10
%One Seg on T5 x10
%Two Seg on T5
%T1 on T5
%T2 on T5
%T3 on T5
%T4 on T5

close all
clear all
clc

%% Anatomy Voxelization:
%Anatomyfilename = 'VoxelDataMultiTarget.mat'; %T1234
%Anatomyfilename = 'VoxelDataTarget1_top.mat'; %T1
%Anatomyfilename = 'VoxelDataTarget2_right.mat'; %T2
%Anatomyfilename = 'VoxelDataTarget3_bottom.mat'; %T3
%Anatomyfilename = 'VoxelDataTarget4_left.mat'; %T4
Anatomyfilename = 'VoxelDataTask2Target1.mat'; %T5
Voxels = load(Anatomyfilename,'Voxel_data');

%% Dexterity analysis settings:
sample_size = 10*1e6; 

%% Set up parallel pool
%Set up the pool until it works
not_worked = true;
while not_worked
    try
        poolobj = parpool('local',[2 30],'SpmdEnabled',false,'IdleTimeout',60);
        disp(poolobj)
        not_worked = false;
        disp('Created Parallel Pool Successfully')
    catch
        disp('Parallel Pool failed retrying...')
        not_worked = true;
        pause(5)
    end
end

%% Rigid Tool
rigid_resultsT5 = zeros(1,10);
for ii=1:10
    %Rigid Tool
    disp(['Starting Fitness function test ' num2str(ii)])
    disp('Testing Design: Rigid Tool on T5')

    tic
    Score_rigid = FastFitnessFunctionRigidTool(sample_size,Anatomyfilename);
    toc

    disp('Dexterity score for This design was:')
    disp(Score_rigid)
    rigid_resultsT5(ii) = Score_rigid;
end
disp('Dexterity scores for all Rigid Tool Analyses:')
disp(rigid_resultsT5)
save('SnakeT5_1_10_results/Rigid_tool_T5_baseline','rigid_resultsT5');

%% One Segment Design 

design.alpha = 0.8700;
design.n = 3;
design.d = 1.1200;
design.w = 4;

directory = 'SnakeT5_1_10_results';

one_resultsT5 = zeros(1,10);
for ii=1:10
    disp(['Starting Fitness function test ' num2str(ii)])
    disp(['Testing a ',num2str(length(design.alpha)),' Segment Design on T5:'])
    disp(design)

    tic
    Score = FastFitnessFunctionVariableSegmentSnakeRobot(design,sample_size,Voxels,directory);
    toc

    disp('Dexterity score for This design was:')
    disp(Score)
    one_resultsT5(ii) = Score;
end
disp('Dexterity scores for Best One Segment Design at T5:')
disp(one_resultsT5)
save('SnakeT5_1_10_results/One_seg_T5_baseline','one_resultsT5');

%% Two Segment Design 

design.alpha = [0.2000 0.8800];
design.n = [3 3];
design.d = [1 1];
design.w = [4 4];

directory = 'SnakeT5_1_10_results';

Two_resultsT5 = zeros(1,10);
for ii=1:10
    disp(['Starting Fitness function test ' num2str(ii)])
    disp(['Testing a ',num2str(length(design.alpha)),' Segment Design on T5'])
    disp(design)

    tic
    Score = FastFitnessFunctionVariableSegmentSnakeRobot(design,sample_size,Voxels,directory);
    toc

    disp('Dexterity score for This design was:')
    disp(Score)
    Two_resultsT5(ii) = Score;
end
disp('Dexterity scores for Best Two Segment Design at T5:')
disp(Two_resultsT5)
save('SnakeT5_1_10_results/Two_seg_T5_baseline','Two_resultsT5');

%% T1 Design 

design.alpha = [0.1300 1.1300];
design.n = [5 3];
design.d = [1.3500 1.0600];
design.w = [4 4];

directory = 'SnakeT5_1_10_results';

T1_resultsT5 = zeros(1,10);
for ii=1:10
    disp(['Starting Fitness function test ' num2str(ii)])
    disp('Testing a T1 Design on T5:')
    disp(design)

    tic
    Score = FastFitnessFunctionVariableSegmentSnakeRobot(design,sample_size,Voxels,directory);
    toc

    disp('Dexterity score for This design was:')
    disp(Score)
    T1_resultsT5(ii) = Score;
end
disp('Dexterity scores for Best T1 Design at T5:')
disp(T1_resultsT5)
save('SnakeT5_1_10_results/T1_on_T5_baseline','T1_resultsT5');

%% T2 Design 

design.alpha = [0.0500 0.6700];
design.n = [6 4];
design.d = [2.8900 1];
design.w = [4 4];

directory = 'SnakeT5_1_10_results';

T2_resultsT5 = zeros(1,10);
for ii=1:10
    disp(['Starting Fitness function test ' num2str(ii)])
    disp('Testing a T2 Design on T5:')
    disp(design)

    tic
    Score = FastFitnessFunctionVariableSegmentSnakeRobot(design,sample_size,Voxels,directory);
    toc

    disp('Dexterity score for This design was:')
    disp(Score)
    T2_resultsT5(ii) = Score;
end
disp('Dexterity scores for Best T2 Design at T5:')
disp(T2_resultsT5)
save('SnakeT5_1_10_results/T2_on_T5_baseline','T2_resultsT5');

%% T3 Design 

design.alpha = [0.2600 1.0100];
design.n = [1 2];
design.d = [7.9400 1.0100];
design.w = [4 4];

directory = 'SnakeT5_1_10_results';

T3_resultsT5 = zeros(1,10);
for ii=1:10
    disp(['Starting Fitness function test ' num2str(ii)])
    disp('Testing a T3 Design on T5:')
    disp(design)

    tic
    Score = FastFitnessFunctionVariableSegmentSnakeRobot(design,sample_size,Voxels,directory);
    toc

    disp('Dexterity score for This design was:')
    disp(Score)
    T3_resultsT5(ii) = Score;
end
disp('Dexterity scores for Best T3 Design at T5:')
disp(T3_resultsT5)
save('SnakeT5_1_10_results/T3_on_T5_baseline','T3_resultsT5');

%% T4 Design 

design.alpha = [0.2900 1.1000];
design.n = [1 3];
design.d = [2.3800 1];
design.w = [4 4];

directory = 'SnakeT5_1_10_results';

T4_resultsT5 = zeros(1,10);
for ii=1:10
    disp(['Starting Fitness function test ' num2str(ii)])
    disp('Testing a T4 Design on T5:')
    disp(design)

    tic
    Score = FastFitnessFunctionVariableSegmentSnakeRobot(design,sample_size,Voxels,directory);
    toc

    disp('Dexterity score for This design was:')
    disp(Score)
    T4_resultsT5(ii) = Score;
end
disp('Dexterity scores for Best T4 Design at T5:')
disp(T4_resultsT5)
save('SnakeT5_1_10_results/T4_on_T5_baseline','T4_resultsT5');

%% T5 Design

% design.alpha = [0.9500 0.9900];
% design.n = [3 2];
% design.d = [9.1500 1.0400];
% design.w = [4 4];

%% END SCRIPT
%End parallel loop delete the pool object
delete(poolobj)
