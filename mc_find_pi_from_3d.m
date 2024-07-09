% This Matlab code performs the Monte-Carlo (MC) simulation to calculate 
% pi (=3.1415926...) number using a concept of geometric probability: 
% probability that point inside sphere is equal to a ratio of volume of a sphere divided by a volume of a cube, 
% where a sphere is inside a cube.      
%
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% July 8, 2024 & University of North Dakota 
%
function [] = mc_find_pi_from_3d
clc;
format long
%
N = 100000;                % number of points in x (y/z) axis. 
%
x = rand(N,1);             % random numbers in [0, 1]
y = rand(N,1);
z = rand(N,1);
%
r2 = (x.^2 + y.^2 + z.^2);   % radius of sphere is one
%
count = 0.;                  % number of point inside the circle 
for i = 1:N
    if (r2(i) <= 1.)
        count = count + 1;
    else
    end
end
pi_estimation = 6 * count/N; % <= relation = volume of sphere/volume of cube

[pi_estimation]

%  N      [pi_estimation]
% 1000    2.976000000000000
% 10000   3.117000000000000
% 100000  3.140460000000000
% 1000000 3.147612000000000

%%%
return
end
