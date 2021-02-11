% Plot the evolutions over time in one graph
% and the boxplots together too...

% Andrew Razjigaev 21 Jan 2021

clc;
clear;
close all;

%% Load data
load('one_seg_fitness'); %Snake1
load('two_seg_fitness'); %Snake2
load('T1_seg_fitness'); %SnakeT1
load('T2_seg_fitness'); %SnakeT2
load('T3_seg_fitness'); %SnakeT3
load('T4_seg_fitness'); %SnakeT4
load('T5_seg_fitness'); %SnakeT5

%% Create Combined Evolutions over time One module to two Module

figure;
%max generations
max_iter = 200;
data = {one_seg_fit, two_seg_fit};
col = {'b','g'};
tex = {'m=1','m=2'};
offset = [-0.0025,0.003];

for ii = length(data):-1:1
    %Extract
    Upper_line = data{ii}.u;
    Lower_line = data{ii}.l;
    mean_fitness = data{ii}.m;
    
    %Extend length of vectors:
    if length(Upper_line)~=max_iter
        Upper_line = [Upper_line; Upper_line(end)*ones(max_iter-length(Upper_line),1)];
        Lower_line = [Lower_line; Lower_line(end)*ones(max_iter-length(Lower_line),1)];
        mean_fitness = [mean_fitness; mean_fitness(end)*ones(max_iter-length(mean_fitness),1)];
    end
    
    %Plot
    generations = linspace(1,length(Upper_line),length(Upper_line))';
    plot(generations,Upper_line,strcat(col{ii},'-'),'LineWidth', 2); hold on;
    plot(generations,Lower_line,strcat(col{ii},'-'),'LineWidth', 2); hold on;
    patch([generations; fliplr(generations')'],[Lower_line; fliplr(Upper_line')'],col{ii},'facealpha',0.1)%[0.7 0.7 1]
    plot(mean_fitness,strcat(col{ii},'--'),'LineWidth', 2); hold on;
    text(max_iter-25,Lower_line(end)+offset(ii),tex{ii},'FontSize',14)
end

title('Mean Fitness and Standard Deviation Over Time');
xlabel('Iteration');
ylabel('Best Dexterity Fitness');
grid on;

%% Create Combined Evolutions over time T1 to T5

figure;
%max generations
max_iter = 250;
data = {T1_seg_fit, T2_seg_fit, T3_seg_fit, T4_seg_fit, T5_seg_fit};
col = {'b','g','r','c','m','y','w'};
tex = {'T1','T2','T3','T4','T5'};
offset = [-0.01,0.015,0.015,-0.01,-0.01];

for ii = length(data):-1:1
    %Extract
    Upper_line = data{ii}.u;
    Lower_line = data{ii}.l;
    mean_fitness = data{ii}.m;
    
    %Extend length of vectors:
    if length(Upper_line)~=max_iter
        Upper_line = [Upper_line; Upper_line(end)*ones(max_iter-length(Upper_line),1)];
        Lower_line = [Lower_line; Lower_line(end)*ones(max_iter-length(Lower_line),1)];
        mean_fitness = [mean_fitness; mean_fitness(end)*ones(max_iter-length(mean_fitness),1)];
    end
    
    %Plot
    generations = linspace(1,length(Upper_line),length(Upper_line))';
    plot(generations,Upper_line,strcat(col{ii},'-'),'LineWidth', 2); hold on;
    plot(generations,Lower_line,strcat(col{ii},'-'),'LineWidth', 2); hold on;
    patch([generations; fliplr(generations')'],[Lower_line; fliplr(Upper_line')'],col{ii},'facealpha',0.1)%[0.7 0.7 1]
    plot(mean_fitness,strcat(col{ii},'--'),'LineWidth', 2); hold on;
    text(max_iter-25,Lower_line(end)+offset(ii),tex{ii},'FontSize',14)
end

title('Mean Fitness and Standard Deviation Over Time');
xlabel('Iteration');
ylabel('Best Dexterity Fitness');
grid on;


%% BOXPLOT Load data

load('one_seg_boxplot'); %Snake1
load('two_seg_boxplot'); %Snake2
load('T1_seg_boxplot'); %SnakeT1
load('T2_seg_boxplot'); %SnakeT2
load('T3_seg_boxplot'); %SnakeT3
load('T4_seg_boxplot'); %SnakeT4
load('T5_seg_boxplot'); %SnakeT5

%% Create Combined Boxplots
alpha_bounds = [0.01 pi/2];
n_bounds = [1 10];
d_bounds = [1 10];

figure;

subplot(2,3,1)
%Create Boxplots ALPHA1
M2 = two_seg_box.alpha1;
T1 = T1_seg_box.alpha1;
T2 = T2_seg_box.alpha1;
T3 = T3_seg_box.alpha1;
T4 = T4_seg_box.alpha1;
T5 = T5_seg_box.alpha1;

group = [ones(size(M2));
        2*ones(size(T1));
        3*ones(size(T2));
        4*ones(size(T3));
        5*ones(size(T4));
        6*ones(size(T5));];
    
   
boxplot([M2; T1; T2; T3; T4; T5],group)
set(gca,'XTickLabel',{'M2','T1','T2','T3','T4','T5'})
ylim(alpha_bounds)
title('Proximal \alpha_1')

subplot(2,3,2)
%Create Boxplots N1
M2 = two_seg_box.n1;
T1 = T1_seg_box.n1;
T2 = T2_seg_box.n1;
T3 = T3_seg_box.n1;
T4 = T4_seg_box.n1;
T5 = T5_seg_box.n1;

group = [ones(size(M2));
        2*ones(size(T1));
        3*ones(size(T2));
        4*ones(size(T3));
        5*ones(size(T4));
        6*ones(size(T5));];
    
   
boxplot([M2; T1; T2; T3; T4; T5],group)
set(gca,'XTickLabel',{'M2','T1','T2','T3','T4','T5'})
ylim(n_bounds)
title('Proximal n_1')

subplot(2,3,3)
%Create Boxplots D1
M2 = two_seg_box.d1;
T1 = T1_seg_box.d1;
T2 = T2_seg_box.d1;
T3 = T3_seg_box.d1;
T4 = T4_seg_box.d1;
T5 = T5_seg_box.d1;

group = [ones(size(M2));
        2*ones(size(T1));
        3*ones(size(T2));
        4*ones(size(T3));
        5*ones(size(T4));
        6*ones(size(T5));];
    
   
boxplot([M2; T1; T2; T3; T4; T5],group)
set(gca,'XTickLabel',{'M2','T1','T2','T3','T4','T5'})
ylim(d_bounds)
title('Proximal d_1')

subplot(2,3,4)
%Create Boxplots ALPHA2
M1 = one_seg_box.alpha2;
M2 = two_seg_box.alpha2;
T1 = T1_seg_box.alpha2;
T2 = T2_seg_box.alpha2;
T3 = T3_seg_box.alpha2;
T4 = T4_seg_box.alpha2;
T5 = T5_seg_box.alpha2;

group = [ones(size(M1));
        2*ones(size(M2));
        3*ones(size(T1));
        4*ones(size(T2));
        5*ones(size(T3));
        6*ones(size(T4));
        7*ones(size(T5));];
    
   
boxplot([M1; M2; T1; T2; T3; T4; T5],group)
set(gca,'XTickLabel',{'M1','M2','T1','T2','T3','T4','T5'})
ylim(alpha_bounds)
title('Distal \alpha_2')

subplot(2,3,5)
%Create Boxplots N2
M1 = one_seg_box.n2;
M2 = two_seg_box.n2;
T1 = T1_seg_box.n2;
T2 = T2_seg_box.n2;
T3 = T3_seg_box.n2;
T4 = T4_seg_box.n2;
T5 = T5_seg_box.n2;

group = [ones(size(M1));
        2*ones(size(M2));
        3*ones(size(T1));
        4*ones(size(T2));
        5*ones(size(T3));
        6*ones(size(T4));
        7*ones(size(T5));];
    
   
boxplot([M1; M2; T1; T2; T3; T4; T5],group)
set(gca,'XTickLabel',{'M1','M2','T1','T2','T3','T4','T5'})
ylim(n_bounds)
title('Distal n_2')

subplot(2,3,6)
%Create Boxplots D2
M1 = one_seg_box.d2;
M2 = two_seg_box.d2;
T1 = T1_seg_box.d2;
T2 = T2_seg_box.d2;
T3 = T3_seg_box.d2;
T4 = T4_seg_box.d2;
T5 = T5_seg_box.d2;

group = [ones(size(M1));
        2*ones(size(M2));
        3*ones(size(T1));
        4*ones(size(T2));
        5*ones(size(T3));
        6*ones(size(T4));
        7*ones(size(T5));];
    
   
boxplot([M1; M2; T1; T2; T3; T4; T5],group)
set(gca,'XTickLabel',{'M1','M2','T1','T2','T3','T4','T5'})
ylim(d_bounds)
title('Distal d_2')