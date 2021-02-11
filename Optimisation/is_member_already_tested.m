function [logic,previous_cost,previous_time] = is_member_already_tested(member_vector,All_populations,All_cost,All_time)
% Given the current population to evolution history, check if the member
% vector has already been computed before:

%Default answer for unique members
logic = false;
previous_cost = [];
previous_time = [];

%Check through Population history
Generations = size(All_populations,2);
population = size(All_populations,1);

%Check each population member in this generation:
for ii = 1:Generations
    for jj = 1:population
        %Check for the case when the member already tested and return its
        %results
        test_member = All_populations{jj,ii};
        if ~isempty(test_member)&&all(member_vector==test_member)
            logic = true;
            previous_cost = All_cost{jj,ii};
            previous_time = All_time{jj,ii};
            break;
        end
    end  
end
end

% %Check through the population history
% Generations = length(All_populations);
% 
% for ii = 1:Generations
%     %Check each population member in this generation:
%     population = length(All_populations(ii).pop);
%     
%     for jj = 1:population
%         %Check for the case when the member already tested and return its
%         %results
%         old_member = All_populations(ii).pop(jj).Position;
%         if ~isempty(old_member)&&all(member_vector==old_member)
%             logic = true;
%             previous_cost = All_populations(ii).pop(jj).Cost;
%             previous_time = All_populations(ii).pop(jj).Time;
%             break;
%         end
%     end  
% end