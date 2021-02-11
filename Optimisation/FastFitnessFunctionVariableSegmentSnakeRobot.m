function [Score] = FastFitnessFunctionVariableSegmentSnakeRobot(design,sample_size,Voxels,directory)
% The fitness is based on the voxelisation of the anatomy and completing a
% dexterity calculation for that. That voxelisation needs to be done in
% GenerateVoxelisation.m first to export the structure of voxel data for
% this.
%author: Andrew Razjigaev 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Joint space limits for design object:%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define the number of segments in this design
segments = length(design.alpha);

%Home Vector:
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
    qipL = -theta_max; qipU = theta_max;
    %tilt
    qitL = -theta_max; qitU = theta_max;
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
V = Voxels.Voxel_data;

%Extract data 
EntranceFrame = V.VoxelBaseFrame;
Ng = V.NumberGoalVoxels;
SSparams = V.ServiceSphere_params;
dV = V.voxelsize;

%Tool Endeffector transfrom to tip
tooltransform = txyz(0,0,5); % 5 mm straight tool on the endeffector

%Determine the Trajectory sampling value based on voxel size and this design:
[traj_length,~] = OptimalTrajLength(EntranceFrame,design,tooltransform,dV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now Complete Dexterity Analysis:%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Reduction variable for parallel loop
ss_map = V.sphere_maps;

parfor ii = 1:N
    %for each configuration q in sampled configuration space Q
    q = Q(ii,:);
    
    %Find Forward kinemtatic endeffector position relative to Entrance
    tend = FastForwardKinematicsSnake(tooltransform,design,q);
    
    %Locate voxel indices
    v = Points2Voxels(V,tend');
    
    %if tend voxel is in bounds (not 0 indices) and labelled goal 
    if (~(any(v==0)) && (V.Goal_labels(v(1),v(2),v(3)) == true))
        
        %Find trajectory, endeffector rotation
        [Traj,Rend,~] = ForwardKinematicsVariableSegmentTraj(EntranceFrame,tooltransform,design,q,dV,traj_length);
        
        %Check if Trajectory is collision free as well
        if (CheckCollision(V,Traj) == false)
          %Get patch in the new map
          new_map = servicesphere_mapping(Rend,V,v);

          %Update sphere maps with OR operation its parfor loop friendly
          ss_map = ss_map|new_map; 
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate total dexterity for the goal%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Average Dexterity = sum(No)/(Ng*Ntheta*Nh);
%Get coordinates in on sphere surface from ss_params
Ntheta = SSparams(1); Nh = SSparams(2); 
%Total dexterity for Design
Total_dexterity = sum(ss_map(:))/(Ng*Ntheta*Nh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Maximum dexterity and distribution %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Determining the dexterity distribution
NoI = subarray_sums(ss_map,[Ntheta,Nh,1]);
dext_dist = NoI./(Ntheta*Nh);

%Maximum Dexterity value and index:
[Max_Dexterity,index] = max(dext_dist(:));
[a,b,c] = ind2sub(size(dext_dist),index);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sum the objectives to get the result %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ouput Total score:
Score = Total_dexterity;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Results as a .mat file for this design %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results = struct('Design_params',design,...
    'mean_dexterity',Total_dexterity,...
    'max_dexterity',Max_Dexterity,...
    'max_location',[a,b,c],...
    'dexterity_distribution',dext_dist,...
    'service_sphere_maps',ss_map);
%Save the Results based on the design parameters
result_file = strcat('Design_alpha',...
    strrep(strrep(strrep(strrep(num2str(design.alpha),'.',''),'        ','_'),'___','_'),'__','_'),...
    '_n',strrep(strrep(num2str(design.n),'.',''),'  ','_'),...
    '_d',strrep(strrep(strrep(strrep(num2str(design.d),'.',''),'        ','_'),'___','_'),'__','_'));

%Just in case something went wrong that is the file already exists 
%while it already exists keep changing the name until its unique
original_result_file = result_file;
duplicate_number = 0;
while exist(strcat(directory,'/',result_file,'.mat'),'file')
    disp('Duplicate Design Name... renaming')
    %This file exists rename it
    duplicate_number = duplicate_number + 1;
    result_file = strcat(original_result_file,'_duplicate_',num2str(duplicate_number));
end

%Ideal way save
%Save in the directory folder
%save(strcat(directory,'/',result_file),'-struct','Results');

%Practical way
%Save in the directory folder
%If failure retry save until it works
not_worked = true;
while not_worked
    try
        %Ideal way save
        save(strcat(directory,'/',result_file),'-struct','Results');
        not_worked = false;
    catch
        disp('Save failed retrying...')
        not_worked = true;
        pause(1) %waits for connection error to resolve
    end
end

end


%{
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
      end
    end


    %Find Forward Kinematic trajectory, endeffector rotation and position
    [Traj,Rend,tend] = ForwardKinematicsVariableSegmentTraj(EntranceFrame,tooltransform,design,q,dV,traj_length);
    
    %Locate voxel indices relative to entrance frame instead of origin
    v = Points2Voxels(V,TransformPoints(Entanceframe2origin,tend'));
    
    %if tend voxel is labelled goal and the shape is collision free 
    if ((V.Goal_labels(v(1),v(2),v(3)) == true) && (CheckCollision(V,Traj) == false))
          %Get patch in the new map
          new_map = servicesphere_mapping(Rend,V,v);

          %Update sphere maps with OR operation its parfor loop friendly
          ss_map = ss_map|new_map; 
    end
%}
