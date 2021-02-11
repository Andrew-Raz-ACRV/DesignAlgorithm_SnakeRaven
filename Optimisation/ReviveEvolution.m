function ReviveEvolution(file)
% Revives a complete evolution file and continues running it for a set 
% amount of generations
%3/8/2020 Andrew Razjigaev

%% Create Directory on HPC-Drive for all the results to go into:

%Add path to all possible folders required
addpath(pwd);
addpath('Snake2_1_10_results');
addpath('SnakeT1_1_10_results');
addpath('SnakeT2_1_10_results');
addpath('SnakeT3_1_10_results');
addpath('SnakeT4_1_10_results');
addpath('SnakeT5_1_10_results');

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

%% Reload the last results:
%file = 'SnakeT5_6.mat';
%for another
iterations = 50; %generations
pre = load(file);

disp('Reviving Evolution file:')
disp(file)
disp(['for ' num2str(iterations) ' generations'])

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

%nVar=3;            % Number of Decision Variables
nVar=6;            % One or Two segment variables 

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

MaxIt=pre.Max_Iterations;      % Maximum Number of Iterations/gnerations

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

Set up the pool until it works
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

%% Revive previous population:
%Reverse Allpop matrix array back to cell array
disp('Starting Snake Evolution Algorithm Reviving Previous Population:');

Allpop = mat2cell(pre.populations_history,ones(1,nPop),nVar*ones(1,MaxIt+1));
Allcost = mat2cell(pre.costs_history,ones(1,nPop),ones(1,MaxIt+1));
Alltime = mat2cell(pre.time_history,ones(1,nPop),ones(1,MaxIt+1));

rng('shuffle'); 

%Recreate Population 100:
empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Time=[];
pop=repmat(empty_individual,nPop,1);

disp('Iteration 100 is regenerating...');
for i=1:nPop 
    % initial population within bound %unifrnd(VarMin,VarMax,VarSize);
    pop(i).Position = Allpop{i,pre.Max_Iterations+1};
    pop(i).Cost = Allcost{i,pre.Max_Iterations+1};
    pop(i).Time = Alltime{i,pre.Max_Iterations+1};
end

%Create Cell array recording all populations over time
Allpop = [Allpop cell(nPop,iterations)];
Allcost = [Allcost cell(nPop,iterations)];
Alltime = [Alltime cell(nPop,iterations)];

%Count the number of function evaluations and repeats
repeat = 0;
func_iter = 0;

disp(['Iteration ' num2str(pre.Max_Iterations) ' has been revived...']);

BestCost=[pre.BestCost; zeros(iterations,1)];
BestSol = pre.BestSol;

%% Continue the Evolution algorithm:

for it=(MaxIt+1):(MaxIt+iterations)
    
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
    
    %Consider Saving data?
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
end
