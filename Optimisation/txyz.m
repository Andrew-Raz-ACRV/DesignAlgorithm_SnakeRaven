function T = txyz(x,y,z)
%Outputs a Transformation matrix for an x, y, z translation
% T = txyz(x,y,z)
%
% That is 
% T = [1 0 0 x;
%      0 1 0 y;
%      0 0 1 z;
%      0 0 0 1];
%
% author: Andrew Razjigaev 2018

T = [eye(3) [x y z]';
     zeros(1,3) 1];
end