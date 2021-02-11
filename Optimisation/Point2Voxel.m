function [ii,jj,kk] = Point2Voxel(Voxel_data,t)
%Given the Voxel space, it finds the closest voxel coordinate for 
%position t in reference to Voxel Base frame, if outside space its 0,0,0.

%Get position in reference to corner of lowest voxel {1,1,1}
%ref = Voxels{1,1,1}.origin;
ref = [Voxel_data.vx(1), Voxel_data.vy(1), Voxel_data.vz(1)];

x = t(1)- ref(1); y = t(2) - ref(2); z = t(3) - ref(3);

%Extract Voxel size
shape = Voxel_data.voxelsize; 
dx = shape(1); dy = shape(2); dz = shape(3);

%Round voxel coordinate to whole voxel i.e. if 0.2 or 0.5 round to 1, 
ii = ceil(x/dx);
jj = ceil(y/dy);
kk = ceil(z/dz);

%Ensure coordinate is within voxelization bounds
bounds = Voxel_data.bounds;

if (ii<1)||(ii>bounds(1))||(jj<1)||(jj>bounds(2))||(kk<1)||(kk>bounds(3))
    ii = 0; jj = 0; kk = 0;
end

end