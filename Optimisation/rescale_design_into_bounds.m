function [p] = rescale_design_into_bounds(y,VarMin,VarMax)
% Saturates the design parameters into its bounds defined in VarMin and
% VarMax
p = y;
for ii = 1:length(y)
    %Saturate
    if p(ii)<VarMin(ii)
        p(ii) = VarMin(ii);
    end
    if y(ii)>VarMax(ii)
        p(ii) = VarMax(ii);
    end
    
    %Round them again:
    
    %round to two decimals
     p(1,ii) = round(p(1,ii),2);
     
     %n must be an integer
     if ii==2||ii==5
         p(1,ii) = round(p(1,ii));
     end
end
end