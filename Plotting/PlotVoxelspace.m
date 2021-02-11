function PlotVoxelspace(Voxel_data)
%Given Voxel_data structure of voxels from Generate Voxelisation
%this plots that voxel space
%Andrew Razjigaev 2019

    dx = Voxel_data.voxelsize(1);
    dy = Voxel_data.voxelsize(2);
    dz = Voxel_data.voxelsize(3);
    vx = Voxel_data.vx;
    vy = Voxel_data.vy;
    vz = Voxel_data.vz;
    EntranceFrame = Voxel_data.VoxelBaseFrame;
    
    %Plot Voxels
    for ii = 1:Voxel_data.bounds(1)
        for jj = 1:Voxel_data.bounds(2)
            for kk = 1:Voxel_data.bounds(3)
                if (Voxel_data.Obstacle_labels(ii,jj,kk)==true)
                    hold on
                    %Plot Voxel  extra logic &&(jj>10)&&(kk<20)
                    plotprism(EntranceFrame,vx(ii),vy(jj),vz(kk),dx,dy,dz,'r')              
                elseif Voxel_data.Goal_labels(ii,jj,kk)==true
                    hold on
                    %Plot Voxel
                    plotprism(EntranceFrame,vx(ii),vy(jj),vz(kk),dx,dy,dz,'g')
                end
            end
        end
    end
    hold on
    plotcoord3(EntranceFrame,10,'r','g','b')
end