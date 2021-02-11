function [Tend] = PlotBasicRigidTool(EntranceFrame,q,dV)
% With the Rigid Tool having the static entrance frame and the joint
% positions q = [pan,tilt,roll,z];
% dV is the voxel size
% This plots the trajectory of the tool

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

%Plot Start of tool
plotcoord3(Tend,2,'r','g','b')
%Draw the insertion tube trajectory
p0 = Tend(1:3,4);
hold on

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
    %Plot Trajectory
    p = Tend(1:3,4);
    plot3([p0(1) p(1)],[p0(2) p(2)],[p0(3) p(3)],'k-');
    plot3(p(1),p(2),p(3),'k*');  
    p0 = p;
    hold on
end

%From the motions: append the tool transforms
%Place n2 points along the z translation
for ii = (n1+1):(n1+n2)
    z_dist = z_offset/n2;
    %Append z translation:
    Tend = Tend*txyz(0,0,z_dist);
    %Append trajectory
    Traj(:,ii) = Tend(1:3,4);
    %Plot Trajectory
    p = Tend(1:3,4);
    plot3([p0(1) p(1)],[p0(2) p(2)],[p0(3) p(3)],'k-');
    plot3(p(1),p(2),p(3),'k*');  
    p0 = p;
    hold on
end

%Curved part sample n3 along curve
dtheta = theta/n3;

for ii = (n1+n2+1):(n1+n2+n3)
    Tend = Tend*Tcurve(r,dtheta);
    %Append trajectory
    Traj(:,ii) = Tend(1:3,4);
    %Plot Trajectory
    p = Tend(1:3,4);
    plot3([p0(1) p(1)],[p0(2) p(2)],[p0(3) p(3)],'k-');
    plot3(p(1),p(2),p(3),'k*');  
    p0 = p;
    hold on
end

%Place n4 points along the z tanslation
for ii = (n1+n2+n3+1):n
    z_dist = z_distal/n4;
    %Append z translation:
    Tend = Tend*txyz(0,0,z_dist);
    %Append trajectory
    Traj(:,ii) = Tend(1:3,4);
    %Plot Trajectory
    p = Tend(1:3,4);
    plot3([p0(1) p(1)],[p0(2) p(2)],[p0(3) p(3)],'k-');
    plot3(p(1),p(2),p(3),'k*');  
    p0 = p;
    hold on
end

%Plot endeffector frame
plotcoord3(Tend,2,'r','g','b')
end
