% This Matlab code performs the Monte-Carlo (MC) simulation and computes a
% following integral using a Metropolis algorithm:
%
%   infinity                 
%    /                    
%    | x^4 * exp(-x^2) dx 
%    /                    
%  -infinity
%  __________________________  = 0.75
%  infinity                    
%    /                    
%    | exp(-x^2) dx 
%    /                    
%  -infinity
%  
%
% Ref. [1] E. Curotto, "Stochastic Simulations of Clusters: Quantum Methods in Flat and Curved Spaces", CPC Press (2010).
%
%      
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% July 8, 2024 & University of North Dakota 
%
function [] = metropolis_example_1
clc; close all; 
format short
%
n_moves = 3000000; % number of Metropolis moves
delta_x = 4.50;   % parameter used to control a rejection rate between 40 and 60%
%
x = 0.0;          % initial coordinate
n_rej = 0.;       % number of rejection 
sm_v = 0.;        % sum of V
sm_v2 = 0.;       % sum of V^2
%
for ii = 1:n_moves
    xt = x + delta_x * (rand(1) - 0.5);
    %
    q = exp(-xt^2)/exp(-x^2); % probability ratio, q = exp(-V(xt))/exp(-V(x))
    if (rand(1) < q)
        x = xt;
    else
        n_rej = n_rej + 1;
    end
    %
    sm_v = sm_v + x^4;
    sm_v2 = sm_v2 + x^4 * x^4;        

end
%
x4_ave = sm_v/n_moves;                            % numerical value
x4_exact = 0.75;                                  % exact value  
x4_sq_ave = sm_v2/n_moves;
sigma_std = sqrt((x4_sq_ave - x4_ave^2)/n_moves); % standard deviation
rejection_ratio = 100*(n_rej/n_moves);            % rejection in percent
%
[x4_ave, x4_exact, sigma_std, rejection_ratio ]
%
%[x4_ave, x4_exact, sigma_std, rejection_ratio ]
% 0.7506    0.7500    0.0014   52.8406


%%%
return
end




