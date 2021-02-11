function PlotBasicOneSegment(EntranceFrame,tooltransform,design,q)
%Plots a line with dots to make a simply representation of a snake like robot
%in 3D plots

%Extract Design values
alpha = design.alpha;
w = design.w;
n = design.n;
d = design.d;

%Compute The Raven Insertion
Traven = [Rx(q(1))*Ry(q(2)) [0 0 0]'; 0 0 0 1]*txyz(0,0,q(3));
Tend = EntranceFrame*Traven;
%Plot the Entrance frame
%plotcoord3(EntranceFrame,2,'r','g','b')
%Plot the Raven frame
%plotcoord3(Tend,2,'r','g','b')
%Draw the insertion tube trajectory
p0 = EntranceFrame(1:3,4);
p = Tend(1:3,4);
plot3([p0(1) p(1)],[p0(2) p(2)],[p0(3) p(3)],'k-');
p0 = p;

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
for ii = 1:n
    if isodd(ii)
        Tend = Tend*Tpan;
        p = Tend(1:3,4);
        plot3([p0(1) p(1)],[p0(2) p(2)],[p0(3) p(3)],'k-');
        plot3(p(1),p(2),p(3),'k*');
    else
        Tend = Tend*Ttilt;
        p = Tend(1:3,4);
        plot3([p0(1) p(1)],[p0(2) p(2)],[p0(3) p(3)],'k-');
        plot3(p(1),p(2),p(3),'k*');
    end
    p0 = p;
end

%Plot Tooltip point
Tend = Tend*tooltransform;
p = Tend(1:3,4);
plot3([p0(1) p(1)],[p0(2) p(2)],[p0(3) p(3)],'k-');
plot3(p(1),p(2),p(3),'k*');
plotcoord3(Tend,2,'r','g','b')

end