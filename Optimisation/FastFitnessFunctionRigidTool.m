function [Score] = FastFitnessFunctionRigidTool(sample_size,Anatomyfilename)
% The fitness is based on the voxelisation of the anatomy and completing a
% dexterity calculation for that. That voxelisation needs to be done in
% GenerateVoxelisation.m first to export the structure of voxel data for
% this. This is for the RIGID TOOL CASE.
%author: Andrew Razjigaev 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Joint space limits for Rigid Tool  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Hand Joint space limits for Rigid Tool:
% %pan x
% qrxL = -pi/6; %radians
% qrxU = pi/6; 
% %tilt y 
% qryL = -pi/6;
% qryU = pi/6;
% %roll z
% qrzL = -pi;
% qrzU = 0; 
% %insertion
% qztL = -30;
% qztU = 20; %mm

%Hand Joint space limits for Rigid Tool:
%pan x
qrxL = -pi; %radians
qrxU = pi; 
%tilt y 
qryL = -pi/2 -pi/4;
qryU = -pi/2 +pi/4;
%roll z
qrzL = -pi/4;
qrzU = +pi/4; 
%insertion
qztL = -50;
qztU = 50; %mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Sampling the Configuration Space  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Joint Space Sampling
N = sample_size;

%Random Joint vectors
qrx = random('unif',qrxL,qrxU,N,1); %N by 1 vector
qry = random('unif',qryL,qryU,N,1);
qrz = random('unif',qrzL,qrzU,N,1);
qzt = random('unif',qztL,qztU,N,1);

%Configuration space calculation where Q(i,:) is the ith configuration in
%Q space
Q = [qrx, qry, qrz, qzt];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now Complete Dexterity Analysis:%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Reduction variable for parallel loop
ss_map = V.sphere_maps;

parfor ii = 1:N
    %for each configuration q in sampled configuration space Q
    q = Q(ii,:);
    %disp(ii);
    
    %Find Forward Kinematic trajectory, endeffector rotation and position
    [Traj,Rend,tend] = ForwardKinematicsRigidToolCurved(EntranceFrame,Entanceframe2origin,q,dV);
    
    if CheckCollision(V,Traj) == false
       %If configuration has no collision or out of bounds in trajectory
       
       %Locate voxel and see if endeffector reached a goal
       t_end = TransformPoints(Entanceframe2origin,tend');
       v = Points2Voxels(V,t_end);

       %Check the endeffector is in bounds (only necessary for the FK here)
       if any(v(:)==0)==false        
          %Check if voxel is goal
          if V.Goal_labels(v(1),v(2),v(3)) == true

              %Get patch in the new map
              new_map = servicesphere_mapping(Rend,V,v);

              %Update sphere maps
              ss_map = ss_map|new_map; 
          end
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
design.alpha = 'RigidTool';
Results = struct('Design_params',design,...
    'mean_dexterity',Total_dexterity,...
    'max_dexterity',Max_Dexterity,...
    'max_location',[a,b,c],...
    'dexterity_distribution',dext_dist,...
    'service_sphere_maps',ss_map);
%Save the Results

%Ensure there are no duplications and ensure that the save is successful:
directory = 'DOF_Rigid_tool_analysis';
result_file = 'Design_RigidTool_analysis';
original_result_file = result_file;

%while it already exists keep changing the name until its unique
duplicate_number = 0;
while exist(strcat(directory,'/',result_file,'.mat'),'file')
    disp('Duplicate Design Name... renaming')
    %This file exists rename it
    duplicate_number = duplicate_number + 1;
    result_file = strcat(original_result_file,'_duplicate_',num2str(duplicate_number));
end

%Save the results:
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