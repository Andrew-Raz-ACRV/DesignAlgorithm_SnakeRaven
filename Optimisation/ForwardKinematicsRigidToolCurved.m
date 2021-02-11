function [Traj,Rend,tend] = ForwardKinematicsRigidToolCurved(EntranceFrame,Entanceframe2origin,q,dV)
% With the Rigid Tool having the static entrance frame and the joint
% positions q = [pan,tilt,roll,z];
% This computes the forward kinematics and trajectory towards the
% endeffector

%Extract the joint motions:
pan = q(1); tilt = q(2); roll = q(3); z = q(4);

%Start at entrance frame:
Tend = EntranceFrame;
%Append x rotation:
Tend = Tend*[Rx(pan) [0 0 0]'; [0 0 0 1]];
%Append y rotation:
Tend = Tend*[Ry(tilt) [0 0 0]'; [0 0 0 1]];
%Append z rotation:
Tend = Tend*[Rz(roll) [0 0 0]'; [0 0 0 1]];


%Design constants:
z_offset = 0;
%Curved part
a = 5; 
b = 20;
r = (b^2 + a^2)/(2*a);
theta = asin(b/r);
%Tool part:
z_distal = 25;

%Sampling trajectory:

%samples along the z translation
n1 = ceil(abs(z)/dV(3));
%samples along the straight part
n2 = ceil(z_offset/dV(3));
%samples along the curved part, the sector arc/voxel
n3 = ceil(theta*r/dV(3));
%samples along the tool
n4 = ceil(z_distal/dV(3));

%Total cascade of data
n = n1 +n2 +n3 +n4;

%Start measuring Trajectory
Traj = zeros(3,n);

%Place n1 points along the z tanslation
for ii = 1:n1
    z_dist = z/n1;
    %Append z translation:
    Tend = Tend*txyz(0,0,z_dist);
    %Append trajectory
    Traj(:,ii) = Tend(1:3,4);
end

%From the motions: append the tool transforms
%Place n2 points along the z translation
for ii = (n1+1):(n1+n2)
    z_dist = z_offset/n2;
    %Append z translation:
    Tend = Tend*txyz(0,0,z_dist);
    %Append trajectory
    Traj(:,ii) = Tend(1:3,4);
end

%Curved part sample n3 along curve
dtheta = theta/n3;

for ii = (n1+n2+1):(n1+n2+n3)
    Tend = Tend*Tcurve(r,dtheta);
    %Append trajectory
    Traj(:,ii) = Tend(1:3,4);
end

%Place n4 points along the z tanslation
for ii = (n1+n2+n3+1):n
    z_dist = z_distal/n4;
    %Append z translation:
    Tend = Tend*txyz(0,0,z_dist);
    %Append trajectory
    Traj(:,ii) = Tend(1:3,4);
end

% Remove the trajectory points that start before the insertion to stay in
% voxel bounds:

% make it relative to entranceframe
Traj_Entrance = (TransformPoints(Entanceframe2origin,Traj'))';

%logic to filter z<0 points from the insertion:
logic = Traj_Entrance(3,:)>0;

%Extract the rows of the trajectory
x_row = Traj(1,:);
y_row = Traj(2,:);
z_row = Traj(3,:);

%Filter each row by the logic
Traj = [x_row(logic); y_row(logic); z_row(logic)];

%End
Rend = Tend(1:3,1:3);
tend = Tend(1:3,4);

end






