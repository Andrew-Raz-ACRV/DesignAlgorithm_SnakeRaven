function [p] = generate_random_design(VarMin,VarMax)
    % generates a gene in bounds VarMin VarMax vectors
    % where VarMin must be equal to size to VarMax
    
    p = zeros(size(VarMax));
    
    for ii = 1:length(VarMax)
        
     %get random value
     p(1,ii) = unifrnd(VarMin(ii),VarMax(ii));
     
     %round to two decimals
     p(1,ii) = round(p(1,ii),2);
     
     %n must be an integer
     if ii==2||ii==5
         p(1,ii) = round(p(1,ii));
     end
    end
    
end