function plotcylinder(rad,d,x0,y0,z0)
%Plots a cylinder of radius, diameter and starting point
%

    %plot cylinder:
    %theta = -2*pi:0.1:2*pi;
    theta = linspace(-2*pi,2*pi,100);
    x = rad*cos(theta)+x0;
    y = rad*sin(theta)+y0;
    zd = z0+d;

    for ii = 1:(length(theta)-1)
        X = [x(ii) x(ii+1) x(ii+1) x(ii)];
        Y = [y(ii) y(ii+1) y(ii+1) y(ii)];
        Z = [z0 z0 zd zd];
        %Z = [zbot(ii) zbot(ii+1) ztop(ii+1) ztop(ii)];
        patch(X,Y,Z,'r');
        hold on
    end

end