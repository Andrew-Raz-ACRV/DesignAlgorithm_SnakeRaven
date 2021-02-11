function PlotDexterityDistribution(Voxel_data,dext_dist,bin)
% Given the Voxel_data structure and the dexterity distribution matrix
% this plots that voxel space with coloured goal voxels corresponding to
% the dexterity value for that voxel. The colour shades depend on 'bin'
% which is a scalar of the number of colour shades to have. e.g. 12
% author: Andrew Razjigaev 2019

%Extract all parameters for voxel plotting
dx = Voxel_data.voxelsize(1);
dy = Voxel_data.voxelsize(2);
dz = Voxel_data.voxelsize(3);
vx = Voxel_data.vx;
vy = Voxel_data.vy;
vz = Voxel_data.vz;
EntranceFrame = Voxel_data.VoxelBaseFrame;

%Colour coding bins: 12 colours
%Create Color from Dexterity
offset = 1;
nc = bin-offset;

%Plot Voxels
for ii = 1:Voxel_data.bounds(1)
    for jj = 1:Voxel_data.bounds(2)
        for kk = 1:Voxel_data.bounds(3)
            if Voxel_data.Goal_labels(ii,jj,kk)==true
                hold on
                Dexterity = dext_dist(ii,jj,kk);
                %Plot Goal Voxel
                Cg = round((nc-1-2*offset)*Dexterity/0.13+1+offset); %magic number??? round(nc*100*Dexterity);
                %Plot Voxel goal with colour code:
                plotprism(EntranceFrame,vx(ii),vy(jj),vz(kk),dx,dy,dz,Cg)             
            end
        end
    end
end

end