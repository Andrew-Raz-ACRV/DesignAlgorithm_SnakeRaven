function PlotBasicVariableSegment(EntranceFrame,tooltransform,design,q)
%Plots a line width representation of a variable segment snake like robot
%in a 3D plot. 
%
% PlotBasicVariableSegment(EntranceFrame,tooltransform,design,q)
%
% Given:
% EntranceFrame = 4x4 Homogeneous Transform that is where the robot starts
% tooltransform = 4x4 Homogeneous Transform that describes the rigid
%                 transform from the last disk to the tip of the end-effector
%
% design = structure containing the segment parameters: alpha, n, d and w
%          To concatenate a two or N segment design define alpha as
%          such: design.alpha = [a1 a2 ... aN]
%          Do the same definition for design.n, design.d and design.w
%
% q = joint vector for the robot defined as:
%     q = [q1 q2 q3 q4 q5 ... qN] joints for N degrees of freedom
%     joints q1,q2,q3 are used to describe the initial transform relative
%     to Entrance frame before defining the snake arm transforms. 
%     For the RAVEN II these are 
%     q1 is x axis rotation Rx(q1)
%     q2 is y axis rotation Ry(q2)
%     q3 is z axis translation txyz(0,0,q3)
%     
%     The other joints are for the segments pan and tilt angles
%     q4, q5 = the first segment pan angle q4 and tilt q5
%     q6, q7 = the second segment pan angle q6 and tilt q7
%
%     In general for k segments the angles are selected following that
%     pattern:
%     pan angle = q(4+(k-1)*2)
%     tilt angle = q(5+(k-1)*2)
%
% author: Andrew Razjigaev 2018

%Extract Design values
alphas = design.alpha; %[a1 a2 a3 a4...]
ws = design.w;
ns = design.n;
ds = design.d;

%Compute The Raven Insertion
%Get RCM Base Frame from Anatomy frame:
isright = true;
RCM_frame = EntranceFrame; %*[Ry(deg2rad(90)) [0 0 0]'; 0 0 0 1];
%Traven = EntranceFrame*[Rx(q(1))*Ry(q(2)) [0 0 0]'; 0 0 0 1]; 

%Compute The Raven Insertion Transform:
La12 = deg2rad(75); La23 = deg2rad(52);

%Select appropriate arm kinematics
if isright==true
    %Right arm:
    T1 = DHmatrix(pi,0,0,q(1));
    T2 = DHmatrix(La12,0,0,q(2));
    T3 = DHmatrix(La23,0,q(3),-pi/2);
else
    %Left arm:
    T1 = DHmatrix(0,0,0,q(1));
    T2 = DHmatrix(La12,0,0,q(2));
    T3 = DHmatrix(pi-La23,0,q(3),pi/2);
end

%Start Kinematic Chain:
Tend = RCM_frame*T1*T2*T3;
% Traven = [Rx(q(1))*Ry(q(2)) [0 0 0]'; 0 0 0 1]*txyz(0,0,q(3));
% baseframe = EntranceFrame*Traven;
% Tend = baseframe;

%Plot the Entrance frame
%plotcoord3(EntranceFrame,2,'r','g','b')
%Plot the Raven frame
%plotcoord3(Tend,2,'r','g','b')
%Draw the insertion tube trajectory
p0 = EntranceFrame(1:3,4);
p = Tend(1:3,4);
plot3([p0(1) p(1)],[p0(2) p(2)],[p0(3) p(3)],'k-','LineWidth',2);
p0 = p;


%Default rules: if this is true then the initial disk is a pan then a tilt
pan_first = true;


%Compute Trajectory while getting forward Kinematics

%Number of segments is length of alpha or w or n or d variables:
segs = length(alphas);
%Initialise the trajectory:
Traj = zeros(3,sum(ns)+1);

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
    
    
    %Plot all disks in a for loop:
    for ii = 1:n
        if pan_first==true %Pan first case
            if isodd(ii) %pan disk for ii = 1,3,5,7...
                Tend = Tend*Tpan;
            else %tilt disk for ii = 2,4,6,8...
                Tend = Tend*Ttilt;
            end           
        else %tilt first case:
            if isodd(ii) %tilt disk for ii = 1,3,5,7...
                Tend = Tend*Ttilt;  
            else %pan disk for ii = 2,4,6,8...
                Tend = Tend*Tpan;
            end    
        end
        %Polt each disk
        p = Tend(1:3,4);
        plot3([p0(1) p(1)],[p0(2) p(2)],[p0(3) p(3)],'k-','LineWidth',2);
        plot3(p(1),p(2),p(3),'k*');  
        hold on
        p0 = p;
    end
    
    %Transition disk addition
    if k~=segs 
        %Embed the inbetween joint with the previous segment
        zeta = deg2rad(-90/segs); %make rotation given cables of equal spacing
        %Transform proximal to distal segment 
        Tpd = [Rz(zeta) [0 0 0]'; 0 0 0 1]; 
        Tend = Tend*Tpd; %Transition transform after nth dik
        %Polt transition disk
        p = Tend(1:3,4);
        plot3([p0(1) p(1)],[p0(2) p(2)],[p0(3) p(3)],'k-','LineWidth',2);
        plot3(p(1),p(2),p(3),'k*');  
        hold on
        p0 = p;
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
end

    %Extract Endeffector final Rotation and position
    %Plot Tooltip point
    Tend = Tend*tooltransform;
    p = Tend(1:3,4);
    plot3([p0(1) p(1)],[p0(2) p(2)],[p0(3) p(3)],'k-','LineWidth',2);
    plot3(p(1),p(2),p(3),'k*');
    plotcoord3(Tend,2,'r','g','b')
    %disp(Tend)
end
