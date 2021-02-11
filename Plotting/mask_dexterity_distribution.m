function [mean_dexterity,dext_dist] = mask_dexterity_distribution(result_filename,Anatomyfilename)
%Given the anatomy Anatomyfilename of Target I, this masks the dexterity distribution in
%result_filename (made for multiple targets) and computes the mean
%dexterity for Target I

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load the Results and Anatomy Data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load Results:
Results = load(result_filename);

Voxels = load(Anatomyfilename);
V = Voxels.Voxel_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Get Mask from Anatomy      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ntheta = V.ServiceSphere_params(1); Nh = V.ServiceSphere_params(2); 
mask = false(size(V.sphere_maps)); %initialise the mask as zeros
%Go through the voxelspace
for ii=1:V.bounds(1)
    for jj=1:V.bounds(2)
        for kk=1:V.bounds(3)
            %If it is a goal select all service sphere patches as true:
            if V.Goal_labels(ii,jj,kk)==true
                a = (ii-1)*Ntheta+1;
                b = (jj-1)*Nh+1;
                mask(a:a+Ntheta-1,b:b+Nh-1,kk) = true(Ntheta,Nh);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Get Dexterity Distribution    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ng = V.NumberGoalVoxels;

ss_map = Results.service_sphere_maps & mask; %mask the map

mean_dexterity = sum(ss_map(:))/(Ng*Ntheta*Nh); % calculate dexterity in new map

%Determining the dexterity distribution
NoI = subarray_sums(ss_map,[Ntheta,Nh,1]);
dext_dist = NoI./(Ntheta*Nh); %segmented dexterity distribution

end