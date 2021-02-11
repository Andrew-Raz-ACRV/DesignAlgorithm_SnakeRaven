function plotDesignSpace()
clc
clf
figure(1)
title('Dexterity Design Space of the Snake Robot')
xlabel('Parameter Alpha Angle')
ylabel('Parameter "n" Number of Joints')
zlabel('Parameter "d" Joint separation distance')
colorbar('Location', 'EastOutside', 'YTickLabel',...
{'0%','1%', '2%','3%','4%','5%','6%','7%','8%','9%','10%','11%','12%','13%','14%','15%'})
view(3)
grid on


alpha = [0.3 1.3100 0.9400 1.3100 1.4300 1]; 
d     = [2 5.1900 8.4100 5.1900 1.4300 1.5]; 
n     = [6 0 9 0 6 7]; 
Output = [-0.0585837 -0.0035 0 -0.0035 -0.0013 -0.1];

% %Design Colour Layers
% ColourScale = linspace(0,10,15);

for ii = 1:length(alpha)
    %Create Color from Dexterity
    nc = 11;
    offset = 1;
    response = -100*Output(ii);
    c = round((nc-1-2*offset)*response/10+1+offset);

    scatter3(alpha(ii), n(ii), d(ii), 30, c,'filled') 
    hold on
    pause(1)
end


title('Dexterity Design Space of the Snake Robot')
xlabel('Parameter Alpha Angle')
ylabel('Parameter "n" Number of Joints')
zlabel('Parameter "d" Joint separation distance')
colorbar('Location', 'EastOutside', 'YTickLabel',...
{'0%','1%', '2%','3%','4%','5%','6%','7%','8%'})

end