function PlotResultsForFitnessFunction(result_filename,Anatomyfilename)
% This function plots all data from the resulting fitness function
% evaluation for a design
% author: Andrew Razjigaev 2019

disp('Visualising the Fitness Function Result Data:')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load the Results and Anatomy Data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load Results:
Results = load(result_filename);

Voxels = load(Anatomyfilename);
V = Voxels.Voxel_data;

%Extract data 
SSparams = V.ServiceSphere_params;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot a rendering of the Result design %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if all(Results.Design_params.alpha(1)=='R')
    %Rigid Tool Case
    
    %Load image of the actual tool for it:
    image = imread('RigidTool_noLogo.png');
    
    %Display the image in a figure:
    figure
    clf
    imshow(imresize(image,1/10))
    title('A Surgical Instrument');
    
else
    %Snake robot Case:
    %Output the Mean dexterity
    disp('SnakeBot Design Parameters:')
    disp(Results.Design_params)

    %lets use a straight configuration:
    segs = length(Results.Design_params.alpha);
    q = zeros(1,3+2*segs);
%     %Straight configuration
%     q(1) = deg2rad(-39.5967);
%     q(2) = deg2rad(-77.9160);

    %Tool Endeffector transfrom to tip
    tooltransform = txyz(0,0,5); % 5 mm straight tool on the endeffector

    %Display a rendering of the snake robot in the next figure 
    figure
    clf
    PlotRenderVariableSegment(eye(4),tooltransform,Results.Design_params,q);
    light('Position',[-1 -1 0.5],'Style','infinite')
    daspect([1 1 1])
    view([-25 15]);
    grid on
    title('Render of the Snake robot Design');
    xlabel('x  - mm');
    ylabel('y  - mm');
    zlabel('z  - mm');
end

%Save Figure
print(gcf,'foo.png','-dpng','-r300');

%Output the Mean dexterity
disp('Mean Dexterity for this design is:')
disp(Results.mean_dexterity)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Maximum Dexterity Service Sphere %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
clf
plotservicesphere(Results.service_sphere_maps,SSparams,Results.max_location)
view([135 -20])
title(['Maximum Dexterity is: ',num2str(100*round(Results.max_dexterity,3)),'%'])         
disp('Max_Dexterity')
disp(Results.max_dexterity)
disp('At Voxel: ')
disp(Results.max_location)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Dexterity Distribution Voxel map %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plotting Dexterity Distribution')
figure
%figure('Name','Task Space','units','normalized','outerposition',[0 0 1 1])
clf
bins = 12;  % number of colour shades
grid on
axis equal
hold on
view([-75 20]); %view([-30 55]); %view([168 76]);
PlotDexterityDistribution(V,Results.dexterity_distribution,bins)
%Results.dexterity_distribution
%Plot
title('Task Space Dexterity Distribution');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');    

%colorbar to the right of the figure:
c = colorbar;
c.Label.String = 'Dexterity as a Percentage %';

%Automatically create the tick labels:
c.Limits = [0 100*round(Results.max_dexterity,3)];
c.Ticks = 100.*round(linspace(0,Results.max_dexterity,bins-1),3);
tick_points = num2cell(c.Ticks);

for ii = 1:length(tick_points)
    tick_points{ii} = num2str(tick_points{ii});
end
c.TickLabels = tick_points;

end