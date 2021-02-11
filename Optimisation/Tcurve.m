function T = Tcurve(r,theta)
%Creates a transform matrix for a curve along the x-z plane
% given r and theta definiing the curve. Good for concentric tube robots

%Find the translations:
z_t = r*sin(theta);
x_t = r - r*cos(theta);

%Get the Rotation that is in the y axis
Rot = Ry(theta);

%Overall Transform:
T = [Rot [x_t; 0; z_t]; [0 0 0 1]];

end