function PlotAnatomy(Voxel_data)
% Plots the anatomy given the Voxel_data structure
% Andrew Razjigaev 2020
    fv = stlread(Voxel_data.filename);
    T = [Rx(deg2rad(33))*Ry(deg2rad(7))*Rz(deg2rad(5)) [-1 5 -42.5]'; 0 0 0 1];
    fv.vertices = TransformPoints(T,fv.vertices); 
    patch(fv,'FaceColor',       [0.8 0.8 0.75], ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
    camlight('headlight');
    material('dull');

end