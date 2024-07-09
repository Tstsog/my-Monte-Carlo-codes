% This Matlab code performs the Monte-Carlo (MC) simulation to calculate 
% pi (=3.1415926...) number using a concept of geometric probability: 
% probability that point inside circle is equal to a ratio of area of a circle divided by area of a square, 
% where circle is inside a square.      
%
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% July 8, 2024 & University of North Dakota 
%
function [] = mc_find_pi_from_2d
clc;
format long
%
N = 5000;           % number of points in x (y) axis.
%
x = rand(N,1);      % random numbers in [0, 1]. 
y = rand(N,1);
r2 = (x.^2 + y.^2); % radius of a circle is one. 
%
count = 0.;         % number of point inside the circle
for i = 1:N
    if (r2(i) <= 1.)
        count = count + 1;
    else
    end
end
pi_estimation = 4 * count/N; % <= relation = area of a circle/area of a square = pi*r^2/4*r^2

[N]
[pi_estimation]

%  N     [pi_estimation]
% 100    3.040000000000000
% 1000   3.088000000000000
% 10000  3.154000000000000
% 50000  3.136880000000000
% 100000 3.132120000000000
% 500000 3.143504000000000

%%%
return
end
