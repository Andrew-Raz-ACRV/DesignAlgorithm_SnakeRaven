function plotDesignSpacefromData()
clc

%Make Figure:
figure(1)
clf


%Extract Results
LoadResults; %load('OptimalDesignSearch11Feb2019.dat')
alpha = OptimalDesignSearch12Feb2019.a;
n = OptimalDesignSearch12Feb2019.n;
d = OptimalDesignSearch12Feb2019.d;
Output = OptimalDesignSearch12Feb2019.v;

%Determine colour based on output:
nc = 11;
offset = 1;
response = -100*Output;
c = response - min(response);
c = round((nc-1-2*offset)*c/max(c)+1+offset);

%Make colourful scatter plot
scatter3(alpha, n, d, 30, c,'filled')
colorbar('Location', 'EastOutside', 'YTickLabel',...
{'0%','1%', '2%','3%','4%','5%','6%','7%','8%','9%','10%','11%','12%','13%','14%','15%'})

%Finishing touches:
title('\fontsize{25} Dexterity Design Space of the Snake Robot')
xlabel('\fontsize{20} Parameter Alpha Angle')
ylabel('\fontsize{20} Parameter "n" Number of Joints')
zlabel('\fontsize{20} Parameter "d" Joint separation distance')
view(3)
axis equal
grid on


end