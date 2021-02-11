function Tendeffector = plotEndeffector(T,L)
% Plots the end effector tool at 4x4 transform T
% This looks like a cylinder with a cone tip
% assuming its that shape L is a z-axis distance that describes the length
% of the end-effector from transform T to the tip of the tool
% e.g.
%     L = tooltransform(3,4);
%     Tendeffector = plotEndeffector(Tend,L);
% 
% The output is the tip coordinate frame matrix Tendeffector
%
% author: Andrew Razjigaev 2018

%Parameters
rad = 1.5; %mm
%L = 7; %mm
%default colour
C = [0.6 0.6 0.6];
%points for patch generation
pts = 100;

%plot cylinder:
theta = linspace(-2*pi,2*pi,pts);
x = rad*cos(theta);
y = rad*sin(theta);

for ii = 1:(length(theta)-1)
    X = [x(ii) x(ii+1) x(ii+1) x(ii)];
    Y = [y(ii) y(ii+1) y(ii+1) y(ii)];
    Z = [0 0 L-rad L-rad];
    %Transfrom coordinates of the patch
    P = TransformPoints(T,[X', Y', Z']);
    X = P(:,1); Y = P(:,2); Z = P(:,3);
    patch(X,Y,Z,C,         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud', ...
         'DiffuseStrength', 1,...
         'BackFaceLighting',    'reverselit', ...
         'AmbientStrength', 0.2);
    hold on
end

Tendeffector = T*txyz(0,0,L);

%Plot end effector
for ii = 1:(length(theta)-1)
    X = [x(ii) x(ii+1) 0 0];
    Y = [y(ii) y(ii+1) 0 0];
    Z = [0 0 rad rad];
    %Transfrom coordinates of the patch
    P = TransformPoints(T*txyz(0,0,L-rad),[X', Y', Z']);
    X = P(:,1); Y = P(:,2); Z = P(:,3);
    patch(X,Y,Z,C,         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud', ...
         'DiffuseStrength', 1,...
         'BackFaceLighting',    'reverselit', ...
         'AmbientStrength', 0.2);
    hold on
end

end