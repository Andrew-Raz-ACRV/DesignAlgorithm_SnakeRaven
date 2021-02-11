function new_map = servicesphere_mapping(Rend,V,Vcoord)
%surface = servicesphere(Rend,surface_prior,ss_params)
% Given the rotation matrix and the prior service sphere and surface sphere
% parameters [Ntheta, Nh, dtheta, dh]
% This updates the service sphere coverage of the surface
   
%Note h is height equal intervals from -1 to 1 based on spherical cap area
% Make grids inclusive like so [theta, theta +dth) [theta+dth, theta+2dth)
% so like in [0, 0.08) ... [2*pi-0.08, 2*pi) -> where 2*pi should be turned
% back to 0 grid thats how to deal with singularity 
% in h its -1 to 1 so 0 to 2: [0, dh) ... but it ens with [2-dh, 2]

%Note:
%Small logical data structure for service spheres:
%sphere_maps = false(nx*Ntheta,ny*Nh,nz);

%Where each Nh by Ntheta service matrix relates to voxel i,j,k
% To extract it would be: 
% sphere_maps(   ((i-1)*Ntheta+1):(i*Ntheta), ...
%                ((j-1)*Nh+1):(j*Nh), ...
%                   k)

%Compute the coordinates on sphere:
z = Rend*[0 0 1]';
%Get unit vector version
z = z/norm(z);
%longitudinal angles
theta = atan2(z(2),z(1));
%latitudinal angles
h = z(3);

%Place point on the sphere surface which is representated as a rectangle
new_map = false(size(V.sphere_maps));

%Get coordinates in on sphere surface from V.ServiceSphere_params
Ntheta = V.ServiceSphere_params(1); Nh = V.ServiceSphere_params(2); 
dtheta = V.ServiceSphere_params(3);  dh = V.ServiceSphere_params(4);

%Find relative coordinate to the top left corner of the matrix which 
%corresponds to point (-pi,1) rather than to the centre (0,0):
th_r = theta - (-pi); h_r = - (h - 1);
%Round the relative angles to the nearest grid:
ii = ceil(th_r/dtheta); %column
jj = ceil(h_r/dh); %row

%Singular point handling 0=2*pi therefore
if (ii == 0)
    ii = Ntheta; %make -pi become +pi position
end

if (jj == 0) %h = -1 goes to 1st grid (-1,-1+dh)
    jj = 1;
end

%Check bounds
if (jj<1)||(jj>Nh)||(ii<1)||(ii>Ntheta)
    %do nothing something is wrong with the calculation
    disp('Failure to update service sphere');
else
    %Everything is okay, update the surface at the voxel:
    a = Vcoord(1); b = Vcoord(2); c = Vcoord(3); 
    %Find corresponding region in sphere maps
    new_map((a-1)*Ntheta + ii, (b-1)*Nh + jj,c) = true; %ii is theta, jj is h
end
    
end