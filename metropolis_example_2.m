% This Matlab code performs the Monte-Carlo (MC) simulation and computes a
% an average energy and heat capacity per particle for an ensemble of N 
% noninteracting one-dimensional harmonic ocsillator (ho) at temperature T, using a Metropolis algorithm.
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
function [] = metropolis_example_2
clc; close all; 
format short
%
k  = 1.50;        % force parameter in harmonic potential   
T_temp = 10.000;  % temperature, (T)
%
n_moves = 200000;  % number of Metropolis moves
delta_x = 15.50;   % parameter used to control a rejection rate between 40 and 60%
%
x = 0.0;          % initial coordinate
n_rej = 0.;       % number of rejection 
sm_v = 0.;        % sum of V
sm_v2 = 0.;       % sum of V^2
%
[V_pot] = pot_ho(k,x); % one-dimensional harmonic oscillator 
%
for ii = 1:n_moves
    xt = x + delta_x * (rand(1) - 0.5);
    %
    [V_pot_t] = pot_ho(k,xt);
    %
    q = exp(-V_pot_t/T_temp)/exp(-V_pot/T_temp); % probability ratio, q = exp(-V(xt)/k_B*T)/exp(-V(x)/k_B*T), 
                                                 % where k_B is Boltzmann constant (k_B = 1 in our case) and T is temperature. 
    if (rand(1) < q)
        x = xt;
        V_pot = V_pot_t;
    else
        n_rej = n_rej + 1;
    end
    %
    sm_v = sm_v + V_pot;
    sm_v2 = sm_v2 + V_pot * V_pot;        

end
%
E_ave = sm_v/n_moves;
E_exact = T_temp/2;                             % exact energy per particle, En = k_B*T/2, where k_B is the Boltzmann constant, in our case, k_B = 1.; 
E2_ave = sm_v2/n_moves;
E2_exact = 3.*(T_temp/2)^2;                     % exact energy per particle, En^2 = 3*(k_B*T/2)^2, where k_B is the Bolztmann constant 
Cv_ave = (E2_ave - E_ave^2)/T_temp^2;           % heat capacity
Cv_exact = 0.5;                                 % exact formula for C_v = 0.5*k_B*T 
sigma_std = sqrt((E2_ave - E_ave^2)/n_moves);   % standart deviation 
rejection_ratio = 100*(n_rej/n_moves);          % rejection in percent

%
[T_temp, E_ave, E_exact, E2_ave, E2_exact, Cv_ave, Cv_exact, sigma_std, rejection_ratio ]

%[T_temp,   E_ave,   E_exact,   E2_ave,  E2_exact,  Cv_ave,  Cv_exact, sigma_std, rejection_ratio ]
% 0.2500    0.1269    0.1250    0.0478    0.0469    0.5065    0.5000    0.0004   51.3240
% 0.5000    0.2525    0.2500    0.1925    0.1875    0.5149    0.5000    0.0008   51.0600
% 1.0000    0.5026    0.5000    0.7583    0.7500    0.5058    0.5000    0.0016   47.5185
% 2.0000    1.0040    1.0000    3.0088    3.0000    0.5002    0.5000    0.0032   48.5800
% 2.5000    1.2354    1.2500    4.5499    4.6875    0.4838    0.5000    0.0039   44.6570
% 3.0000    1.5042    1.5000    6.8208    6.7500    0.5065    0.5000    0.0048   50.6355
% 5.0000    2.5061    2.5000   18.6898   18.7500    0.4964    0.5000    0.0079   45.6980
% 8.0000    4.0167    4.0000   48.9847   48.0000    0.5133    0.5000    0.0128   47.0515
%10.0000    5.0004    5.0000   74.7343   75.0000    0.4973    0.5000    0.0158   50.7465



%%%
return
end

%
function [V_pot ] = pot_ho(k,x)
%
V_pot = 0.5*k*x.^2;
%%%
return
end


