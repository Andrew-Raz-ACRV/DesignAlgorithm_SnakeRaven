function PlotRenderOneSegment(EntranceFrame,tooltransform,design,q)
%Plots a rendering of a snake like robot
%in 3D plots

%Extract Design values
alpha = design.alpha;
w = design.w;
n = design.n;
d = design.d;

%Compute The Raven Insertion
Traven = [Rx(q(1))*Ry(q(2)) [0 0 0]'; 0 0 0 1]*txyz(0,0,q(3));
baseframe = EntranceFrame*Traven;
%Plot the Entrance Frame
plotcoord3(baseframe,2,'r','g','b')

%Compute Trajectory while getting forward Kinematics

%Additional parameters
rad = w/2;
r = rad / sin(alpha);
np = round(n/2);
nt = n - np;

%Extract pan and tilt
theta_p = q(4);
theta_t = q(5);

%Pan trasform  0 to 1
r01 = [0 2*r*(1 - cos(alpha)*cos(theta_p/(2*np)))*sin(theta_p/(2*np)) ...
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

%Plot the snake Robot
plotcoord3(baseframe,w,'r','g','b')
Tsnake = baseframe;

Ttiltrot = [Rz(pi/2) [0 0 0]';
            zeros(1,3) 1];

for ii = 1:n
    if isodd(ii)
        Tsnake = Tsnake*Tpan;
        plotVariableNeutralLineDisk(Tsnake,w,alpha,d)
    else
        Tsnake = Tsnake*Ttilt;
        plotVariableNeutralLineDisk(Tsnake*Ttiltrot,w,alpha,d)
    end
    hold on;
end

plotcoord3(Tsnake,w,'r','g','b')
hold on
plotInsertionTube(baseframe,w,alpha,30);
hold on
L = tooltransform(3,4);
Tendeffector = plotEndeffector(Tsnake,L);
hold on
plotcoord3(Tendeffector,w,'r','g','b');
hold on
end