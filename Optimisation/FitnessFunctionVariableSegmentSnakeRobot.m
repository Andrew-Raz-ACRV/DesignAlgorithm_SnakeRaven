function [Score] = FitnessFunctionVariableSegmentSnakeRobot(design,sample_size,plotting,Anatomyfilename)
% The fitness is based on the voxelisation of the anatomy and completing a
% dexterity calculation for that. That voxelisation needs to be done in
% GenerateVoxelisation.m first to export the structure of voxel data for
% this.
%author: Andrew Razjigaev 2019

%Define the number of segments in this design
segments = length(design.alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Joint space limits for design object:%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Raven Joints: Lower Upper
% %x rotation
% qrxL = 0; %radians
% qrxU = pi/4; %pi/4
% %y rotation
% qryL = -pi/4;
% qryU = pi/4;
% %z translation
% qrzL = 0;
% qrzU = 50; %mm

q_home = [deg2rad(-39.5967) deg2rad(-77.9160)]; % right home position
%q_home = [deg2rad(39.5967) deg2rad(-102.084)]; % left home position

%Raven Joints: Lower Upper
%x rotation
qrxL = q_home(1) - pi/4; %radians was 0
qrxU = q_home(1) + pi/4; %was pi/4
%y rotation
qryL = q_home(2) - pi/4;
qryU = q_home(2) + pi/4;
%z translation
qrzL = 0; %50
qrzU = 100; %mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Configuration space:%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Joint Space Sampling
N = sample_size;

%Raven Joint space:
qrx = random('unif',qrxL,qrxU,N,1); %N by 1 vector
qry = random('unif',qryL,qryU,N,1);
qrz = random('unif',qrzL,qrzU,N,1);

%Configuration space calculation where Q(i,:) is the ith configuration in
%Q space
Q = zeros(N,3+segments*2);
Q(:,1:3) = [qrx, qry, qrz];

%Append Configurations for each Segment Joints:
for ii = 1:segments
    %Continuum Joints: Lower Upper
    theta_max = (design.alpha(ii)*design.n(ii))/2; %Maximum bending
    %pan
    qipL = -theta_max;
    qipU = theta_max;
    %tilt
    qitL = -theta_max;
    qitU = theta_max;
    %Sample configuration space, pan and tilt joints:
    qip = random('unif',qipL,qipU,N,1);
    qit = random('unif',qitL,qitU,N,1);
    %Append segment configurations to whole configuration space:
    Q(:,(4+(ii-1)*2):(5+(ii-1)*2)) = [qip, qit];
    % pan = 4+(k-1)*2  :seg =1 4+0, 5+0
    %tilt = 5+(k-1)*2   :seg = 2 4+2 = 6, 5+2 = 7
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load the Voxel Data about Task Space %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Voxels = load(Anatomyfilename);
V = Voxels.Voxel_data;

%Extract data 
EntranceFrame = V.VoxelBaseFrame;
Entanceframe2origin = V.Baseframe2origin;
Ng = V.NumberGoalVoxels;
SSparams = V.ServiceSphere_params;
dV = V.voxelsize;

%Tool Endeffector transfrom to tip
tooltransform = txyz(0,0,5); % 5 mm straight tool on the endeffector

%Determine the Trajectory sampling value based on voxel size and this design:
[traj_length,~] = OptimalTrajLength(EntranceFrame,design,tooltransform,dV);

%Plot Anatomy
if plotting==true
    disp('Showing Workspace...')
    figure('Name','Task Space','units','normalized','outerposition',[0 0 1 1])
    clf
    %Load Anatomy
%     PlotAnatomy(V)
%     hold on
    %Plot Voxels
    disp('Making voxels')
    PlotVoxelspace(V)  
    %Figure settings
    grid on
    axis equal
    title('Task Space');
    view([-30 60]);    %[180 15]); %145 0 shows entrance angle port
    xlabel('X axis');
    ylabel('Y axis');
    zlabel('Z axis'); 
    axis([-15 15 -20 20 -80 10]);
    disp('Made voxels')
    pause(0.001)
    %Display a rendering of the snake robot in the next figure 
    figure(2)
    clf
    q = zeros(size(Q(1,:)));
    PlotRenderVariableSegment(EntranceFrame,tooltransform,design,q)
    light('Position',[-1 -1 0.5],'Style','infinite')
    daspect([1 1 1])
    grid on
    title('Render of the Initial Design of the Snake robot');
    xlabel('x  - mm');
    ylabel('y  - mm');
    zlabel('z  - mm');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now Complete Dexterity Analysis:%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Reduction variable for parallel loop
ss_map = V.sphere_maps;

for ii = 1:N
    %for each configuration q in sampled configuration space Q
    q = Q(ii,:);
    disp(ii);
    
    %Old errorous version:
    %[Traj,Rend,tend] = ForwardKinematicsVariableSegment(EntranceFrame,tooltransform,design,q);
    
    %Find Forward Kinematic trajectory, endeffector rotation and position
    [Traj,Rend,tend] = ForwardKinematicsVariableSegmentTraj(EntranceFrame,tooltransform,design,q,dV,traj_length);
    
    %Plot Configuration onto the voxel space
%     if plotting==true
%         figure(1)
%         %clf
%         %PlotVoxelspace(V) 
%         hold on
%         %Add new stick-like figure of the snake robot
%         PlotBasicVariableSegment(EntranceFrame,tooltransform,design,q)
%         view([-30 60]);
%         grid on
%         axis equal
%         axis([-15 15 -20 20 -80 10]);
%     end

    if CheckCollision(V,Traj) == false
       %If configuration has no collision or out of bounds in trajectory
       
       %Locate voxel and see if endeffector reached a goal
       t_end = TransformPoints(Entanceframe2origin,tend');
       v = Points2Voxels(V,t_end);

       %Check if voxel is goal
      if V.Goal_labels(v(1),v(2),v(3)) == true
          
          %Get patch in the new map
          new_map = servicesphere_mapping(Rend,V,v);

          %Update sphere maps
          ss_map = ss_map|new_map; 
          
          %Plot the service sphere update:
          if plotting==true
               %Plot Successful reach in anatomy
               figure(1)
               hold on
               %Add new stick-like figure of the snake robot
               PlotBasicVariableSegment(EntranceFrame,tooltransform,design,q)
               grid on
               axis equal
               view([-30 60]);
               axis([-15 15 -20 20 -80 10]);
               disp([num2str(ii),' Target was reached successfully!'])
               pause(0.001)
               
%                %plot Update of Sphere
%                figure(3)
%                clf
%                plotservicesphere(ss_map,SSparams,v)
          end
      else
%           if plotting==true
%               disp([num2str(ii),' No collision but no target'])
%           end
      end
    else
%         if plotting==true
%             disp([num2str(ii),' Collision out of bounds occured'])
%         end
    end
end

%Reduction variable for parallel loop
V.sphere_maps = ss_map;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate total dexterity for the goal%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Average Dexterity = sum(No)/(Ng*Ntheta*Nh);
%Get coordinates in on sphere surface from ss_params
Ntheta = SSparams(1); Nh = SSparams(2); 
%Total dexterity for Design
f1 = sum(ss_map(:))/(Ng*Ntheta*Nh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Maximum dexterity and distribution %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Determining the dexterity distribution
NoI = subarray_sums(ss_map,[Ntheta,Nh,1]);
dext_dist = NoI./(Ntheta*Nh);

%Maximum Dexterity value and index:
[Max_Dexterity,index] = max(dext_dist(:));
[a,b,c] = ind2sub(size(dext_dist),index);

% %Maximum Dexterity surface sphere:
% Max_surface = ss_map((a-1)*Ntheta+1 : a*Ntheta, (b-1)*Nh+1 : b*Nh, c);

%Plot Maximum Dexterity Service Sphere
figure
clf
plotservicesphere(ss_map,SSparams,[a,b,c])
title(['Maximum Dexterity is: ',num2str(Max_Dexterity)])         

disp('Max_Dexterity')
disp(Max_Dexterity)
disp('At Voxel: ')
disp([a,b,c])

%{
%Show max dexterity voxel in anatomy:
show_max = false;

if show_max == true
    disp('Showing Workspace...')
    figure('Name','Task Space','units','normalized','outerposition',[0 0 1 1])
    clf
    %Load Anatomy
%     PlotAnatomy(VoxelizationAnatomy)
%     hold on
    %Goal prism Volume origin and dimensions relative to entrance frame
%     g_ori = [-15,-30,20]';
%     g_dim = [20, 20, 35]';
%     plotprism(EntranceFrame,g_ori(1),g_ori(2),g_ori(3),g_dim(1),g_dim(2),g_dim(3),'g')
    %Plot Voxels
    disp('Making voxels')
    PlotVoxelspace(Voxel_data)  
    hold on
    p = Voxels{a,b,c}.origin;
    plotprism(EntranceFrame,p(1),p(2),p(3),VoxelizationAnatomy.voxelsize(1),VoxelizationAnatomy.voxelsize(2),VoxelizationAnatomy.voxelsize(3),'b')
    hold on
    %Figure settings
    grid on
    axis equal
    title('Task Space');
    view([250 55]);    %[180 15]); %145 0 shows entrance angle port
    xlabel('X axis');
    ylabel('Y axis');
    zlabel('Z axis');     
    pause(1)
    %End plotting
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the simplicity penalty %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha_constraints = [0 pi/2];
% d_constraints = [0 50];
% n_constraints = [1 100];
% 
% %Todo:
% alpha_penalty = exponential_penalty(design.alpha,sigma,alpha_constraints);
% d_penalty = exponential_penalty(design.d,sigma,d_constraints);
% n_penalty = exponential_penalty(design.n,sigma,n_constraints);
% 
% %Total Penalty for chosen design parameters
% f2 = alpha_penalty+d_penalty+n_penalty;
f2 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sum the objectives to get the result %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ouput Total score:
Score = f1 + f2;

if plotting==true
    disp('Dexterity score for This design was:')
    disp(Score)
end

end
