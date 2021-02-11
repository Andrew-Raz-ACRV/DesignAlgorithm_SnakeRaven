function [logic] = is_in_bounds(p,pmin,pmax)
    %calculates whether a vector is inside the minimum amd maximum limits
    % is_in_bounds([0 3 6],[1 1 1],[5 5 5]) = false
    logic = and(all(p<=pmax),all(p>=pmin));
end