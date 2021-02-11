% Differential Evolution script
% SnakeRaven Optimisation on HPC
% Andrew Razjigaev 4 May 2020

clc;
clear;
close all;

%% Create Directory on HPC-Drive for all the results to go into:
addpath(pwd);
%Create a directory name based on the current time
directory = strcat('Snake_Evolution_Results',strrep(strrep(datestr(datetime),':','_'),' ','_'));

%Create new directory ensure it doesn't already exist
[~, msg, ~] = mkdir(directory);
pause(5)

%If it already exists i.e. a message about it keeps appearing, try again until it works
%this happens if the job was run at the same time in seconds as another evolution job
while(~isempty(msg))
    directory = strcat('Snake_Evolution_Results',strrep(datestr(datetime),':','_'));
    [~, msg, ~] = mkdir(directory);
    pause(3)
end

disp('All fitness results in this logfile are being saved in folder:');
disp(directory)

%Add path to directory
addpath(directory);

%% Anatomy Voxelization:
Anatomyfilename = 'VoxelDataMultiTarget.mat';
% Anatomyfilename = 'VoxelDataTarget1_top.mat';
% Anatomyfilename = 'VoxelDataTarget2_right.mat';
% Anatomyfilename = 'VoxelDataTarget3_bottom.mat';
% Anatomyfilename = 'VoxelDataTarget4_left.mat';
%sample_size = 8*1e6; 
sample_size = 10*1e6;  % must be larger than 9100*(18*9) orientation patches
% must be larger than 1474200

disp('Using Anatomy file:')
disp(Anatomyfilename)
disp('Configuration sample size:')
disp(sample_size)

%% Problem Definition

%Voxel Map
Voxels = load(Anatomyfilename,'Voxel_data');

% Cost Function
CostFunction=@(design) FastFitnessFunctionVariableSegmentSnakeRobot(design,sample_size,Voxels,directory);

nVar=3;            % Number of Decision Variables
%nVar=6;            % One or Two segment variables 

VarSize=[1 nVar];   % Decision Variables Matrix Size

% alpha n d bounds: [lower upper]
alpha_bounds = [0.01 pi/2];
n_bounds = [1 10];
d_bounds = [1 10];
% resolution of alpha, n and d
res = [0.01 1 0.01];
if nVar==3
    % One segment
    VarMin=[alpha_bounds(1) n_bounds(1) d_bounds(1)];          % Lower Bound of Decision Variables
    VarMax=[alpha_bounds(2) n_bounds(2) d_bounds(2)];          % Upper Bound of Decision Variables
    disp('Solving a one segment design 3 variables')
    %Solve the design space:
    disp('That makes a design space of this many designs:')
    N_alpha = round((alpha_bounds(2) - round(alpha_bounds(1),2))/res(1)) + 1; % that is 157
    N_n = round((n_bounds(2) - n_bounds(1))/res(2)) + 1; %that is 10
    N_d = round((d_bounds(2) - d_bounds(1))/res(3)) + 1; % that is 901
    Design_space = N_alpha*N_n*N_d;
    disp(Design_space)
elseif nVar==6
    % Two segment
    VarMin=[alpha_bounds(1) n_bounds(1) d_bounds(1) ...
        alpha_bounds(1) n_bounds(1) d_bounds(1)];          % Lower Bound of Decision Variables
    VarMax=[alpha_bounds(2) n_bounds(2) d_bounds(2) ...
        alpha_bounds(2) n_bounds(2) d_bounds(2)];          % Upper Bound of Decision Variables  
    disp('Solving a two segment design 6 variables')
    %Solve the design space:
    disp('That makes a design space of this many designs:')
    N_alpha = round((alpha_bounds(2) - round(alpha_bounds(1),2))/res(1)) + 1; % that is 157
    N_n = round((n_bounds(2) - n_bounds(1))/res(2)) + 1; %that is 10
    N_d = round((d_bounds(2) - d_bounds(1))/res(3)) + 1; % that is 901
    Design_space = N_alpha*N_n*N_d*N_alpha*N_n*N_d;
    disp(Design_space)
end

%% DE Parameters

MaxIt=100;      % Maximum Number of Iterations/gnerations

nPop=10*nVar;        % Population Size

F = 0.5;        % Amplification Factor 0 - 2

pCR=0.8;        % Crossover Probability

disp(['Running for ' num2str(MaxIt) ' iterations/generations'])
disp(['population size: ' num2str(nPop)])
%disp(['Mutation factor range: ' num2str(beta_min) ' to ' num2str(beta_max)])
disp(['Amplification/mutation factor: ' num2str(F)])
disp(['Crossover Probability: ' num2str(pCR)])

%% Set up parallel pool

poolobj = parpool('local',[2 30],'SpmdEnabled',false,'IdleTimeout',60);
disp(poolobj)

%% Initialization

disp('Starting Snake Evolution Algorithm Creating Initial Population:');

rng('shuffle'); 
% Avoids Repeating the same number generator for 
% every instance of this script on HPC
% https://au.mathworks.com/help/matlab/math/why-do-random-numbers-repeat-after-startup.html
% https://au.mathworks.com/help/matlab/ref/rng.html
%
% Initializes generator based on the current time, resulting in a different 
% sequence of random numbers after each call to rng. As jobs are forced to
% have unique folders based on time, each instance will therefore have a
% different sequence of random umbers

empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Time=[];

BestSol.Cost=inf;

pop=repmat(empty_individual,nPop,1);

%Create Cell array recording all populations over time
Allpop = cell(nPop,MaxIt+1);
Allcost = cell(nPop,MaxIt+1);
Alltime = cell(nPop,MaxIt+1);

%Count the number of function evaluations and repeats
repeat = 0;
func_iter = 0;

disp('Iteration 0 has started...');
for i=1:nPop 
    % initial population within bound %unifrnd(VarMin,VarMax,VarSize);
    pop(i).Position=generate_random_design(VarMin,VarMax);
    
    % Run the Fitness Function
    disp(['Testing member ' num2str(i) ' of generation 0']);
    disp('Evaluating design: ')
    disp(vector2designstruct(pop(i).Position))
    
    
    %Check if the design has already been tested:
    [was_tested, prior_cost, prior_time] = is_member_already_tested(pop(i).Position,Allpop,Allcost,Alltime);
    
    if was_tested
        repeat = repeat + 1;
        disp('Already Evaluated skipping recalculation')
        pop(i).Cost = prior_cost;
        pop(i).Time = prior_time;
    else
        %Unique Design needs To be calculated
        func_iter = func_iter + 1;
        tic
        pop(i).Cost=-1*CostFunction(vector2designstruct(pop(i).Position));
        toc
        pop(i).Time = toc;
    end
    
    disp('Dexterity score for This design was:')
    disp(-1*pop(i).Cost)
    
    if pop(i).Cost<BestSol.Cost
        BestSol=pop(i);
    end
    
    %Record population member and data:
    Allpop{i,1} = pop(i).Position;
    Allcost{i,1} = pop(i).Cost;
    Alltime{i,1} = pop(i).Time;
    
end
disp('Iteration 0 has finished...');

BestCost=zeros(MaxIt,1);

%% DE Main Loop

for it=1:MaxIt
    
    disp(['Iteration ' num2str(it) ' has started...']);
    
    for i=1:nPop
        
        x=pop(i).Position; % Get Population member gene
        
        A=randperm(nPop); % Get a random ordering of the population
        
        A(A==i)=[]; %Ensure the order doesn't include the current member
        
        a=A(1);
        b=A(2);
        c=A(3);
        
        % Mutation
        %V = Xr1 + F (Xr2 - Xr3)
        v = pop(a).Position + F.*(pop(b).Position - pop(c).Position);
        
        % rescale into bounds i.e. saturate min and max after mutation
        v = rescale_design_into_bounds(v,VarMin,VarMax);
		
        % Crossover
        u=zeros(size(x));
        j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=pCR
                u(j)=v(j);
            else
                u(j)=x(j);
            end
        end
        
        % rescale into bounds i.e. saturate min and max after crossover
        NewSol.Position=rescale_design_into_bounds(u,VarMin,VarMax);
        
        % Run the Fitness Function
        disp(['Testing member ' num2str(i) ' of generation ' num2str(it)]);
        disp('Evaluating design: ')
        disp(vector2designstruct(NewSol.Position))

        %Check if the design has already been tested:
        [was_tested, prior_cost, prior_time] = is_member_already_tested(NewSol.Position,Allpop,Allcost,Alltime);
    
        if was_tested
            repeat = repeat + 1;
            disp('Already Evaluated skipping recalculation')
            NewSol.Cost = prior_cost;
            NewSol.Time = prior_time;
        else
            %Unique Design needs To be calculated
            func_iter = func_iter + 1;
            tic
            NewSol.Cost=-1*CostFunction(vector2designstruct(NewSol.Position));
            toc
            NewSol.Time = toc;
        end           
    
        disp('Dexterity score for This design was:')
        disp(-1*NewSol.Cost)
        
        %Record New population member and data:
        Allpop{i,it+1} = NewSol.Position;
        Allcost{i,it+1} = NewSol.Cost;
        Alltime{i,it+1} = NewSol.Time;
        
        %Survival of the fittest:
        if NewSol.Cost<pop(i).Cost
            pop(i)=NewSol;
            
            if pop(i).Cost<BestSol.Cost
               BestSol=pop(i);
            end
        end
       
    end
    
    % Update Best Cost
    BestCost(it)=BestSol.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ' finished: Best Cost = ' num2str(-BestCost(it))]);
    disp('\n');
end

%% Show Results

%End parallel loop delete the pool object
delete(poolobj)

%Closing Message:
disp('Evolution time complete. The Best solution was:');
disp(vector2designstruct(BestSol.Position))
disp('With best Dexterity:')
disp(-1*BestSol.Cost)
disp('\n');

%Show some statistics:
disp('In a design space of this many designs:')
disp(Design_space)
disp('Total Fitness evaluations called for evolution:')
disp(nPop * (MaxIt+1))
disp('Number of actual unique designs evaluated:')
disp(func_iter)
disp('Number of actual repeated designs:')
disp(repeat)

%Create Optimal results file:
OptimResults = struct('BestSol',BestSol,...
    'BestCost',BestCost,...
    'Max_Iterations',MaxIt,...
    'populations_history',cell2mat(Allpop),...
    'costs_history',cell2mat(Allcost),...
    'time_history',cell2mat(Alltime));
%Reverse Allpop matrix array back to cell array
%mat2cell(OptimResults.populations_history,ones(1,nPop),nVar*ones(1,MaxIt+1))

%Save the Results
Evolution_file = strcat(directory,' Finished_ ',strrep(strrep(datestr(datetime),':','_'),' ','_'));

%Save Evolutionresults until works:
%Stupid Directory access failure retry save until works
not_worked = true;
cd(directory);
while not_worked
    try
        %save(strcat(directory,'/',Evolution_file),'-struct','OptimResults');
        save(Evolution_file,'-struct','OptimResults');
        not_worked = false;
        disp('Save succssful end of evolution')
    catch
        disp('Save failed retrying...')
        not_worked = true;
        pause(5)
    end
end

