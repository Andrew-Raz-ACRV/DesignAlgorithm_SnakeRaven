function [Traj,Rend,tend] = ForwardKinematicsVariableSegment(EntranceFrame,tooltransform,design,q)
%Computes forward Kinematics for a variable segment continuum on raven insertion
%tube. tend = [x y z]', Rend = rotation matrix. 
%The design alpha = [a1 a2 a3... aN] that is N segment robot
% Traj = [P1, P2, ... tend] where P(i) is each position of disk in robot

%Extract Design values
alphas = design.alpha; %[a1 a2 a3 a4...]
ws = design.w;
ns = design.n;
ds = design.d;

%Compute The Raven Insertion
nr = 10;
dz = q(3)/nr;
Traven = EntranceFrame*[Rx(q(1))*Ry(q(2)) [0 0 0]'; 0 0 0 1];
TrajRaven = zeros(3,nr);
for ii = 1:nr
Traven = Traven*txyz(0,0,dz);
TrajRaven(:,ii) = Traven(1:3,4); 
end
Tend = Traven;

%Old method
% Traven = [Rx(q(1))*Ry(q(2)) [0 0 0]'; 0 0 0 1]*txyz(0,0,q(3));
% baseframe = EntranceFrame*Traven;
% Tend = baseframe;


%Default rules: if this is true then the initial disk is a pan then a tilt
pan_first = true;


%Compute Trajectory while getting forward Kinematics

%Number of segments is length of alpha or w or n or d variables:
segs = length(alphas);
%Initialise the trajectory:
Traj = zeros(3,sum(ns)+1); %3*sum(ns)+5 -new  
%Variables for disk counting
first_disk = 1;

for k = 1:segs
    %Get the design parameters for the first segment:
    alpha = alphas(k);
    w = ws(k);
    n = ns(k);
    d = ds(k);
    
    %Compute other paramters:
    rad = w/2;
    r = rad / sin(alpha);
    
    %Decide how many pan and tilt disks are needed
    if pan_first == true
        %if pan first then split n into np and nt like so       
        np = round(n/2);
        nt = n - np;        
    else
        %if tilt first then split n into np and nt like so       
        nt = round(n/2);
        np = n - nt; 
    end
    
    %Extract Joint angles pan and tilt for kth segment
    theta_p = q(4+(k-1)*2); %i.e k =1 q4, k=2 its q6
    theta_t = q(5+(k-1)*2);
    
    %Compute Transform Matrices
    %Pan trasform  0 to 1
    r01 = [0 -2*r*(1 - cos(alpha)*cos(theta_p/(2*np)))*sin(theta_p/(2*np)) ...
        2*r*(1 - cos(alpha)*cos(theta_p/(2*np)))*cos(theta_p/(2*np))]';
    r11 = [0 0 d]';

    Tpan = [Rx(theta_p/np) r01+Rx(theta_p/np)*r11;
            zeros(1,3)                         1];

    %Tilt  transform 1 to 2
    r12 = [2*r*(1 - cos(alpha)*cos(theta_t/(2*nt)))*sin(theta_t/(2*nt)) 0 ...
        2*r*(1 - cos(alpha)*cos(theta_t/(2*nt)))*cos(theta_t/(2*nt))]';
    r22 = [0 0 d]';

    Ttilt = [Ry(theta_t/nt) r12+Ry(theta_t/nt)*r22;
            zeros(1,3)                         1];
    
    
    %Plot all disks in the segment in a for loop based on total disk length:
    for ii = first_disk:(first_disk+n-1)
        if pan_first==true %Pan first case
            if isodd(ii-(first_disk-1)) %pan disk for ii = 1,3,5,7... relative to first disk
                Tend = Tend*Tpan;
            else %tilt disk for ii = 2,4,6,8...
                Tend = Tend*Ttilt;
            end           
        else %tilt first case:
            if isodd(ii-(first_disk-1)) %tilt disk for ii = 1,3,5,7... relative to first disk
                Tend = Tend*Ttilt;  
            else %pan disk for ii = 2,4,6,8...
                Tend = Tend*Tpan;
            end    
        end
        %Discrete Trajectory is the position of each disk
        Traj(:,ii) = Tend(1:3,4);
        %Traj(:,3*ii) = Tend(1:3,4);
        %Add a few points in between transition
        %Traj(:,3*ii-1) = Tend(1:3,4) - [0 0 -d/1.5]';
        %Traj(:,3*ii-2) = Tend(1:3,4) - [0 0 -d/1.5]';
    end
    
    
    %Transition disk addition is embedded into the nth disk of the segment
    if k~=segs 
        %Embed the inbetween joint with the previous segment
        zeta = deg2rad(-90/segs); %make rotation given cables of equal spacing
        %Transform proximal to distal segment 
        Tpd = [Rz(zeta) [0 0 0]'; 0 0 0 1]; 
        Tend = Tend*Tpd; %Transition transform after nth dik
    end
      
    %Decide if the next segment should be pan first or not based on the pan
    %tilt pattern being preserved.
    if (pan_first==true)&&(isodd(n))
        %e.g. n=3 : p t p -t
        pan_first = false;
    elseif (pan_first==false)&&(isodd(n))
        %e.g. n=3 : t p t -p
        pan_first = true;
    elseif (pan_first==true)&&(isodd(n)==false)
        %e.g. n=4 : p t p t -p
        pan_first = true;
    elseif (pan_first==false)&&(isodd(n)==false)
        %e.g. n=4 : t p t p -t
        pan_first = false;     
    end
    
    %Update first_disk for next segment
    first_disk = first_disk + n;
end

    %Extract Endeffector final Rotation and position
    Tend = Tend*tooltransform;
    Rend = Tend(1:3,1:3);
    tend = Tend(1:3,4);

    %Append tool tip to trajectory remember the tool is 5mm
    Traj(:,end) = tend;
%     Traj(:,end-1) = tend - [0 0 -1]';
%     Traj(:,end-2) = tend - [0 0 -2]';
%     Traj(:,end-3) = tend - [0 0 -3]';
%     Traj(:,end-4) = tend - [0 0 -4]';
%     
%     %Combine Raven and trajectory snake
%     Traj = [TrajRaven, Traj];
end