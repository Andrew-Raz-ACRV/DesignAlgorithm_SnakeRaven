function [logic] = CheckCollision(Voxel_data,points)
% Given the Voxel details, check if points [x,y,z]'
% intersect an obstacle voxel true if collision false for no collision
% all vectorised code

%Default assume it has no collision
logic = false;

%Assume points are relative to global frame make it in reference to voxel
%frame before using NearestVoxel
Points = TransformPoints(Voxel_data.Baseframe2origin,points');

%Find corresponding Voxels: Points = [x, y, z] -> Voxel = [ii, jj, kk]
[Voxelcoords] = Points2Voxels(Voxel_data,Points);

if any(Voxelcoords(:)==0)
    %Its outside the Voxel space bounds
    logic = true;
    %disp('Outside voxel bounds')
    %Assume if outside voxel bounds its a 'collision'
else
    %Voxelcoords is all in the boundaries
    %Check for obstacles: for loop 
    Extracted = false(size(Voxelcoords,1),1);
    for ii = 1:size(Voxelcoords,1)
        Extracted(ii,1) = Voxel_data.Obstacle_labels(Voxelcoords(ii,1),Voxelcoords(ii,2),Voxelcoords(ii,3));
    end
    if any(Extracted)
        logic = true;
        %disp('Collision voxel occured one of these is true')
    end
end
end