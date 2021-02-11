function [combos,number] = combinations(v,n)
% finds all combinations for a code like so:
% you have a scalar 'n' length code with code values defined in vector 'v'
% This function gives you all combinations of the code and the number of
% combinations that exist.
% That means you'll get a matrix of combinations of size numel(v)^n by n
% It also gives the number of combinations that is numel(v)^n
%
% Example Problem 1:
% Find all combinations in a problem of getting values [1 2 3]
% if you pick a value twice:
% i.e. v = [1 2 3], n = 2
% The number of combinations are 3^2 = 9 they are intuitively:
% [1 1, 1 2, 1 3, 2 1, 2 2, 2 3, 3 1, 3 2, 3 3] 
%
% What about picking the number 3 times:
% i.e. n = 3
% 3^3 = 27 combinations
% [1 1 1, 1 1 2, 1 1 3, 1 2 1, 1 2 2, 1 2 3, 1 3 1, ... , 3 3 2, 3 3 3]
%
% Example problem 2:
% you want to brute force a 4 digit pin code for a pinpad that gives 10
% possible digits - 0 to 9 - for a pin code. How many combinations and what
% are those combinations?
%
% Take v as 0:9 and n as 4 to use the function
% [combos, number] = combinations(0:9,4)
%
% number = 10^4 = 10000
% combos =
% [0 0 0 0, 0 0 0 1, 0 0 0 2, 0 0 0 3, ... , 9 9 9 7, 9 9 9 8, 9 9 9 9]
%
  
%Number of combinations
number = numel(v)^n; %Same as number = size(combos,1);

%Use binomial nchoosek to get all the combos
step1 = kron(ones(1,n),v);
step2 = nchoosek(step1,n);
combos = unique(step2,'rows');
   
end