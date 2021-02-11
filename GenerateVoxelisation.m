%Generate Anatomy Voxelisation:
%Andrew Razjigaev 15 June 2020
%Generates the voxel maps for all task spaces
%Shows how difficult those targets are using the obstacle occlusion index

clc
clear all
close all

%% Load Anatomy
tic
%Load Anatomy
Afile = 'kneemodel.stl';
fv = stlread(Afile);

T = [Rx(deg2rad(33))*Ry(deg2rad(7))*Rz(deg2rad(5)) [-1 5 -42.5]'; 0 0 0 1];

fv.vertices = TransformPoints(T,fv.vertices); 

figure('Name','Task Space','units','normalized','outerposition',[0 0 1 1])

% patch(fv,'FaceColor',       [0.8 0.8 0.75], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
% camlight('headlight');
% material('dull');

grid on
axis('image');
%view([-45 45]);
view([-30 55]); %Best viewing angle of all four targets

xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');

hold on 

%% User defined Regions:

%Entrance frame RCM is at Origin! plotcoord3(eye(4),10,'r','g','b')
Entranceframe = [Ry(deg2rad(-90)) [0 0 0]'; 0 0 0 1];
plotcoord3(Entranceframe,10,'r','g','b')

%Labelling the Goal Volume

%MultiGoal:
g_ori_dim =  [-41,-5,-9, 7, 7, 7; %top target
              -60,-2,-7, 6, 6, 6; %Bottom target
              -52, 5,-9, 9, 8, 6; %Left target
              -49,-17,-7,9, 5, 7; %right target
              -48, 12, 1, 10, 7, 10]; % Task 2

%Specific Goals
%  g_ori_dim =  [-41,-5,-9, 7, 7, 7]; %1 top target
%  g_ori_dim =  [-60,-2,-7, 6, 6, 6]; %3 Bottom target
%  g_ori_dim =  [-52, 5,-9, 9, 8, 6]; %4 Left target
%  g_ori_dim =  [-49,-17,-7,9, 5, 7]; %2 right target

%Task 2 Goals
%MultiGoal:
%g_ori_dim =  [-60,-2,-7, 6, 10, 6]; %[-52, 15, 5, 5, 5, 5];
%g_ori_dim =  [-48, 12, 1, 10, 7, 10]; %50.29 official


%Making the bounds
NG = size(g_ori_dim,1);
G_xbounds = zeros(NG,2); 
G_ybounds = zeros(NG,2); 
G_zbounds = zeros(NG,2);
for ii = 1:NG
    G_xbounds(ii,:) = [g_ori_dim(ii,1)  g_ori_dim(ii,1)+g_ori_dim(ii,4)];
    G_ybounds(ii,:) = [g_ori_dim(ii,2)  g_ori_dim(ii,2)+g_ori_dim(ii,5)];
    G_zbounds(ii,:) = [g_ori_dim(ii,3)  g_ori_dim(ii,3)+g_ori_dim(ii,6)];
end

% %Plot Goal Region
% for ii = 1:NG
%     hold on
%     plotprism(Entranceframe,g_ori_dim(ii,1),g_ori_dim(ii,2),g_ori_dim(ii,3),g_ori_dim(ii,4),g_ori_dim(ii,5),g_ori_dim(ii,6),'g')  
% end

%Trocar space
t_ori = [-25, -7, -5]';
t_dim = [25, 12, 12]';
t_xbounds = [t_ori(1) t_ori(1)+t_dim(1)];
t_ybounds = [t_ori(2) t_ori(2)+t_dim(2)];
t_zbounds = [t_ori(3) t_ori(3)+t_dim(3)];

% %Plot Trocar Region
% hold on
% plotprism(Entranceframe,t_ori(1),t_ori(2),t_ori(3),t_dim(1),t_dim(2),t_dim(3),'b')

%% Define Voxelisation Region

%Lets say a voxel space of dxmm*dymm*dzmm cubes space
dx = 2; %was 1 1 1
dy = 2;
dz = 2;

%The size of the voxel space:
L = 30; %25; % z_upper-z_lower
W = 40; % x_upper-x_lower
H = 70; % y_upper-y_lower

%Define scope of Voxel space
%Real space bounds at rounded limits
x_lower = -H; %-round(W/2);
x_upper = 0; %round(W/2);
y_lower = -round(W/2); %-round(H/2)-15;
y_upper = round(W/2); %round(H/2)-15;
z_lower = -round(L/2); %0;
z_upper = round(L/2); %round(L);

%Define Service Sphere parameters:
Ntheta = 72/4; %number of longitudinal segments 
Nh = 36/4; %number of latitudinal segments 
dtheta = 2*pi/Ntheta;
dh = 2/Nh;

%Points of Voxel corners including the final point:
vx = linspace(x_lower,x_upper,round((x_upper-x_lower)/dx)+1)'; 
vy = linspace(y_lower,y_upper,round((y_upper-y_lower)/dy)+1)'; 
vz = linspace(z_lower,z_upper,round((z_upper-z_lower)/dz)+1)';

%Number of Voxels in each axis:
nx = length(vx)-1;
ny = length(vy)-1;
nz = length(vz)-1;

%New Matrix Voxel matrix 

%Initialise labels as logical falses
null_matrix = false(nx,ny,nz);
goal_label = null_matrix;
obst_label = null_matrix;

%i.e. Voxel(i,j,k) real coord centroid is just
% vx(i)+dx/2  vy(j)+dy/2  vz(k)+dz/2
vxcoord = vx(1:end-1)+(dx/2);
vycoord = vy(1:end-1)+(dy/2);
vzcoord = vz(1:end-1)+(dz/2);

%Small logical data structure for service spheres:
sphere_maps = false(nx*Ntheta,ny*Nh,nz);

%Where each Nh by Ntheta service matrix relates to voxel i,j,k
% To extract it would be: 
% sphere_maps(   ((i-1)*Ntheta+1):(i*Ntheta), ...
%                ((j-1)*Nh+1):(j*Nh), ...
%                   k)

% %Plot Voxel Mesh to visualise voxel space
% hold on
% %Plot Grid X-Z lines:
% for jj = 1:length(vy)
%     
%     for ii = 1:length(vz)
%         p = [vx,vy(jj)*ones(size(vx)),vz(ii)*ones(size(vx))];
%         P = TransformPoints(Entranceframe,p); 
%         plot3(P(:,1),P(:,2),P(:,3),'k-','LineWidth',0.1) %plots xz lines
%         hold on
%     end
% 
%     for ii = 1:length(vx)
%         p = [vx(ii)*ones(size(vz)), vy(jj)*ones(size(vz)), vz];
%         P = TransformPoints(Entranceframe,p);
%         plot3(P(:,1),P(:,2),P(:,3),'k-','LineWidth',0.1) %plots zx lines
%         hold on
%     end
% end
% 
% %Plot Grid Y-X lines:
% for jj = 1:length(vz)
%     for ii = 1:length(vx)
%         p = [vx(ii)*ones(size(vy)),vy,vz(jj)*ones(size(vy))];
%         P = TransformPoints(Entranceframe,p);
%         plot3(P(:,1),P(:,2),P(:,3),'k-','LineWidth',0.1) %plots zx lines
%         hold on
%     end
% end

%% Defining the output Structure: Voxel_data

%All Voxel data fields for a Voxel structure
field1 = 'VoxelBaseFrame'; value1 = Entranceframe;
field2 = 'voxelsize'; value2 = [dx dy dz];
field3 = 'XBounds'; value3 = [x_lower x_upper];
field4 = 'YBounds'; value4 = [y_lower y_upper];
field5 = 'ZBounds'; value5 = [z_lower z_upper];
field6 = 'TotalNumberofVoxels'; value6 = nx*ny*nz;
field7 = 'NumberGoalVoxels'; value7 = 0;
field8 = 'NumberObstacleVoxels'; value8 = 0;
field9 = 'NumberEmptyVoxels'; value9 = value6 - value7 - value8;
field10 = 'ServiceSphere_params'; value10 = [Ntheta, Nh, dtheta, dh];
field11 = 'Baseframe2origin'; value11 = inv(Entranceframe);
field12 = 'filename'; value12 = Afile;
field13 = 'vxcoord'; value13 = vxcoord;
field14 = 'vycoord'; value14 = vycoord;
field15 = 'vzcoord'; value15 = vzcoord;
field16 = 'bounds'; value16 = [nx, ny, nz];
field17 = 'Goal_labels'; value17 = goal_label;
field18 = 'Obstacle_labels'; value18 = obst_label;
field19 = 'sphere_maps'; value19 = sphere_maps;
field20 = 'vx'; value20 = vx;
field21 = 'vy'; value21 = vy;
field22 = 'vz'; value22 = vz;
field23 = 'Obstacle_occlusion_index'; value23 = zeros(nx,ny,nz);
field24 = 'OOI_shortestpath'; value24 = zeros(nx,ny,nz);
field25 = 'Combined_OOI'; value25 = zeros(nx,ny,nz);

%Define Voxelisation Structure
Voxel_data = struct(field1,value1,field2,value2,field3,value3,...
    field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,...
    field9,value9,field10,value10,field11,value11,field12,value12,...
    field13,value13,field14,value14,field15,value15,field16,value16,...
    field17,value17,field18,value18,field19,value19,field20,value20,...
    field21,value21,field22,value22,field23,value23,field24,value24,field25,value25);

%% Labelling the Voxels as Obstacles or Goals using tissue location

%Go through each obstacle point and label the nearest voxel as 'obst'

%Obstacle points
obstacle_x = fv.vertices(:,1);
obstacle_y = fv.vertices(:,2);
obstacle_z = fv.vertices(:,3);

% %Filter obstacle points to remove patella voxels
% obstacle_x = obstacle_x(obstacle_z<-30);
% obstacle_y = obstacle_y(obstacle_z<-30);
% obstacle_z = obstacle_z(obstacle_z<-30);

Entanceframe2origin = inv(Entranceframe);

%Labelling of Voxels
for ii = 1:length(obstacle_x)
    %Define Obstacle location:
    obst = [obstacle_x(ii); obstacle_y(ii); obstacle_z(ii)];
    %Make the obstacle location in reference to the entrance frame
    Obst = TransformPoints(Entanceframe2origin,obst');
    
    %Find nearest Voxel
    [a,b,c] = Point2Voxel(Voxel_data,Obst);
    %Check that its a valid coordinate
    if any(isnan([a,b,c]))
        %if its outside neglect it
    else
        %i.e. none of them are nan its inside voxel space
        %Check if the tissue is in Goal volume bounds
        if isinBounds(Obst,G_xbounds,G_ybounds,G_zbounds)
            %Make the voxel be labelled as goal
            Voxel_data.Goal_labels(a,b,c) = true;
        else
            %if its not part of the trocar then its an obstacle:
            if ~isinBounds(Obst,t_xbounds,t_ybounds,t_zbounds)
                Voxel_data.Obstacle_labels(a,b,c) = true;
            end
        end
    end
end

%% Dilate Obstacles by width of the Snake robot
%This is necessary so that we can assume the snake robot has zero thickness

%assume conservative width / diameter of snake robot
w = 4; %mm

%Define voxel side cubic Structuring element of voxels
se_x = ceil(w/dx); %e.g. 2 voxel grid = 4/2
se_y = ceil(w/dy);
se_z = ceil(w/dz);

%Make the structuring element an odd window. i.e. if se = 2 make it 3
if isodd(se_x)==false
    se_x = se_x + 1;
elseif isodd(se_y)==false
    se_y = se_y + 1;
elseif isodd(se_z)==false
    se_z = se_z + 1;
end

%Dilation with logical cubic structuring element:
se = true(se_x,se_y,se_z);

%Dilate 3D voxels
Voxel_data.Obstacle_labels = imdilate(Voxel_data.Obstacle_labels,se);

%Ensure no both 'goal and obstacle' true cases solve with XOR operation:
exclusive_labels = xor(Voxel_data.Obstacle_labels,Voxel_data.Goal_labels);

%and exclusive and goal gives exclusively only goal labels that are not
%obstacles
Voxel_data.Goal_labels = and(Voxel_data.Goal_labels,exclusive_labels);

%% Plot and count Obstacle and Goal Voxels

%Plot Voxels
for ii = 1:nx
    for jj = 1:ny
        for kk = 1:nz
            if Voxel_data.Obstacle_labels(ii,jj,kk)==true
                hold on
                %Plot Obstacle Voxel
                plotprism(Entranceframe,vx(ii),vy(jj),vz(kk),dx,dy,dz,'r')              
            elseif Voxel_data.Goal_labels(ii,jj,kk)==true
                hold on
                %Plot Goal Voxel
                plotprism(Entranceframe,vx(ii),vy(jj),vz(kk),dx,dy,dz,'g')
            end
        end
    end
end

%Counting the labels:
Voxel_data.NumberGoalVoxels = sum(Voxel_data.Goal_labels(:));
Voxel_data.NumberObstacleVoxels = sum(Voxel_data.Obstacle_labels(:));
Voxel_data.NumberEmptyVoxels = Voxel_data.TotalNumberofVoxels - Voxel_data.NumberGoalVoxels - Voxel_data.NumberObstacleVoxels;

%% Obstacle occlusion index Calculation Local

%This calculation tells us how difficult it is to reach the targets
%based on the number of obstacles around each goal voxel
[combos, number] = combinations(-2:2,3); % Get all neighbor points relative to origin
for ii = 1:nx
    for jj = 1:ny
        for kk = 1:nz
            if Voxel_data.Goal_labels(ii,jj,kk)==true
                %Calculate the index there are 26 neighbor voxels
                %Go through all voxel positions in the 3x3x3 cube
                pos = repmat([ii,jj,kk],number,1) + combos; %make it relative to current voxel
                obstacle_sum = 0;
                for s = 1:number
                    a = pos(s,1); b = pos(s,2); c = pos(s,3); 
                    if ~is_in_bounds([a,b,c],[1 1 1],[nx ny nz]) || Voxel_data.Obstacle_labels(a,b,c)
                        % Accumulate the obstacle occlusions
                        obstacle_sum = obstacle_sum + 1;
                    end
                end
                Voxel_data.Obstacle_occlusion_index(ii,jj,kk) = -obstacle_sum/(number-1);
            end
        end
    end
end

%Plot the Obstacle Occlusion Index in a figure
figure;
offset = 1;
bins = 12; %colour shades
nc = bins-offset;
for ii = 1:nx
    for jj = 1:ny
        for kk = 1:nz             
            if Voxel_data.Goal_labels(ii,jj,kk)==true
                hold on
                %calculate colour based on index
                ooi = Voxel_data.Obstacle_occlusion_index(ii,jj,kk);
                Cg = round((nc-1-2*offset)*ooi/0.13+1+offset); %magic number??? round(nc*100*Dexterity);
                %Plot Goal Voxel
                plotprism(Entranceframe,vx(ii),vy(jj),vz(kk),dx,dy,dz,Cg)
            end
        end
    end
end
grid on
axis('image');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
view([-75 20]);
title('Obstacle Occlusion Index: Local Obstacle Count')

%colorbar to the right of the figure:
c = colorbar();
c.Label.String = 'Obstacle Occlusion Index as a Percentage %';

max_value = min(Voxel_data.Obstacle_occlusion_index(:));
%min_value = max(Voxel_data.Obstacle_occlusion_index(:));

%Automatically create the tick labels:
c.Limits = [100*round(max_value,3) 0];
c.Ticks = flip(100.*round(linspace(0,max_value,bins-1),3));
tick_points = num2cell(c.Ticks);

for ii = 1:length(tick_points)
    tick_points{ii} = num2str(tick_points{ii});
end
c.TickLabels = tick_points;

%% Obstacle occlusion index Calculation Global

%Determine number of obstacles in a straight line between entrance frame and goal
%Entrance Point
Pa = Entranceframe(1:3,4); %Entrance point at Origin:
for ii = 1:nx
    for jj = 1:ny
        for kk = 1:nz
            if Voxel_data.Goal_labels(ii,jj,kk)==true
                Pb = [vxcoord(ii); vycoord(jj); vzcoord(kk)]; %Centroid of the voxel
                %Get linear Points From Pa to Pb
                Np = ceil(norm(Pb-Pa)/dx); %number of points
                px = linspace(Pa(1),Pb(1),Np)';
                py = linspace(Pa(2),Pb(2),Np)';
                pz = linspace(Pa(3),Pb(3),Np)';
                %Convert Them to voxel coordinates:
                [Vcoords] = Points2Voxels(Voxel_data,[px, py, pz]);
                Nv = size(Vcoords,1);
                %Count the number of obstacles in those points:
                value = 0;
                for s=1:Nv
                    a = Vcoords(s,1); b = Vcoords(s,2); c = Vcoords(s,3); 
                    if ~is_in_bounds([a,b,c],[1 1 1],[nx ny nz]) || Voxel_data.Obstacle_labels(a,b,c)
                        value = value + 1;
                    end
                end
                Voxel_data.OOI_shortestpath(ii,jj,kk) = -value/Nv; %Percentage of obstacles in the shorest path
            end
        end
    end
end

%Plot the OOI shortest path in a figure
figure;
offset = 1;
bins = 12; %colour shades
nc = bins-offset;
for ii = 1:nx
    for jj = 1:ny
        for kk = 1:nz             
            if Voxel_data.Goal_labels(ii,jj,kk)==true
                hold on
                %calculate colour based on index
                ooi = Voxel_data.OOI_shortestpath(ii,jj,kk);
                Cg = round((nc-1-2*offset)*ooi/0.13+1+offset); %magic number??? round(nc*100*Dexterity);
                %Plot Goal Voxel
                plotprism(Entranceframe,vx(ii),vy(jj),vz(kk),dx,dy,dz,Cg)
            end
        end
    end
end
grid on
axis('image');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
view([-75 20]);
title('Obstacle Occlusion Index: Obstacles in Shortest path')

%colorbar to the right of the figure:
c = colorbar();
c.Label.String = 'Obstacle Occlusion Index as a Percentage %';

max_value = min(Voxel_data.OOI_shortestpath(:));
%min_value = max(Voxel_data.OOI_shortestpath(:));

%Automatically create the tick labels:
c.Limits = [100*round(max_value,3) 0];
c.Ticks = flip(100.*round(linspace(0,max_value,bins-1),3));
tick_points = num2cell(c.Ticks);

for ii = 1:length(tick_points)
    tick_points{ii} = num2str(tick_points{ii});
end
c.TickLabels = tick_points;

%% Combine the Obstacle Occlusion Indices:

%Get the mean of the two results
Voxel_data.Combined_OOI = (Voxel_data.OOI_shortestpath+Voxel_data.Obstacle_occlusion_index)/2;

%Plot the Combined_OOI in a figure
figure;
offset = 1;
bins = 12; %colour shades
nc = bins-offset;
for ii = 1:nx
    for jj = 1:ny
        for kk = 1:nz             
            if Voxel_data.Goal_labels(ii,jj,kk)==true
                hold on
                %calculate colour based on index
                ooi = Voxel_data.Combined_OOI(ii,jj,kk);
                Cg = round((nc-1-2*offset)*ooi/0.13+1+offset); %magic number??? round(nc*100*Dexterity);
                %Plot Goal Voxel
                plotprism(Entranceframe,vx(ii),vy(jj),vz(kk),dx,dy,dz,Cg)
            end
        end
    end
end
grid on
axis('image');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
view([-75 20]);
title('Combined Obstacle Occlusion Index')

%colorbar to the right of the figure:
c = colorbar();
c.Label.String = 'Obstacle Occlusion Index as a Percentage %';

max_value = min(Voxel_data.Combined_OOI(:));
%min_value = max(Voxel_data.Combined_OOI(:));

%Automatically create the tick labels:
c.Limits = [100*round(max_value,3) 0];
c.Ticks = flip(100.*round(linspace(0,max_value,bins-1),3));
tick_points = num2cell(c.Ticks);

for ii = 1:length(tick_points)
    tick_points{ii} = num2str(tick_points{ii});
end
c.TickLabels = tick_points;

%% End of Generate Voxelisation

%Display Structure Data:
disp(Voxel_data)
whos Voxel_data

%Tidy up:
clearvars -except Voxel_data 

%Save the Voxel_data
%get rid of of hyphen in the date
save(strcat('VoxelData',strrep(date,'-','')),'-v7.3')

%savefig(strcat('VoxelData',strrep(date,'-','')))
%saveas(gcf,strcat('VoxelData',strrep(date,'-',''),'png'))
toc