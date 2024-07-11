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
T_temp = 10.0000;  % temperature, (T)
%
n_moves = 500000;  % number of Metropolis moves
delta_x = 8.500;   % parameter used to control a rejection rate between 40 and 60%
%
x = 0.0;          % initial coordinate
n_rej = 0.;       % number of rejection 
sm_v = 0.;        % sum of V
sm_v2 = 0.;       % sum of V^2
%
[V_pot] = pot_ho(k,x); % one-dimensional harmonic oscillator 
%
%%%%%%%%%%%%%%%%%%%%%%%%
fileID_save_data_1 = fopen('metropolis_example_2.txt','w');
%
pot_save = zeros(n_moves, 1);
pot_sq_save = zeros(n_moves, 1);
%
for ii = 1:n_moves
    %
    xt = x + delta_x * (rand(1) - 0.5);
    %
    [V_pot_t] = pot_ho(k,xt);
    %
    q = exp(-V_pot_t/T_temp)/exp(-V_pot/T_temp); % probability ratio, q = exp(-V(xt)/k_B*T)/exp(-V(x)/k_B*T), 
                                                 % where k_B is Boltzmann constant (k_B = 1 in our case) and T is temperature. 
    if (rand(1) < q)
        x = xt;
        V_pot = V_pot_t;
        %
        sm_v = V_pot;
        sm_v2 = V_pot * V_pot;               
    else
        n_rej = n_rej + 1;
    end

%%%
    pot_save(ii) = sm_v;
    pot_sq_save(ii) = sm_v2;
    output = [ii, pot_save(ii), pot_sq_save(ii)];
    %
    fprintf(fileID_save_data_1, '%4.6f \t %4.6f \t %8.6f\n', output); 
end
fclose(fileID_save_data_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_mc_data = fopen('metropolis_example_2.txt', 'r');               % 
read_mc_data = textscan(read_mc_data, '%f %f %f');
mc_step_ii = read_mc_data{1};
mc_pot_val = read_mc_data{2};
mc_pot_sq_val = read_mc_data{3};


%%%
E_ave = sum(mc_pot_val)/length(mc_pot_val);
E_exact = T_temp/2;                             % exact energy per particle, En = k_B*T/2, where k_B is the Boltzmann constant, in our case, k_B = 1.; 
E2_ave = sum(mc_pot_sq_val)/length(mc_pot_sq_val);
E2_exact = 3.*(T_temp/2)^2;                     % exact energy per particle, En^2 = 3*(k_B*T/2)^2, where k_B is the Bolztmann constant 
Cv_ave = (E2_ave - E_ave^2)/T_temp^2;           % heat capacity
Cv_exact = 0.5;                                 % exact formula for C_v = 0.5*k_B*T 
sigma_std = sqrt((E2_ave - E_ave^2)/n_moves);   % standart deviation 
rejection_ratio = 100*(n_rej/n_moves);          % rejection in percent
%
%%%
[T_temp, E_ave, E_exact, E2_ave, E2_exact, Cv_ave, Cv_exact, sigma_std, rejection_ratio, delta_x ]

%%%
%[T_temp,   E_ave,   E_exact,   E2_ave,  E2_exact,  Cv_ave,  Cv_exact, sigma_std, rejection_ratio, delta_x ]
% 0.2500    0.1255    0.1250    0.0471    0.0469    0.5022    0.5000    0.0003   51.4356    2.5000
% 0.5000    0.2496    0.2500    0.1865    0.1875    0.4968    0.5000    0.0005   45.4830    3.0000
% 1.0000    0.5028    0.5000    0.7576    0.7500    0.5048    0.5000    0.0010   39.0448    3.5000
% 5.0000    2.4867    2.5000   18.5241   18.7500    0.4936    0.5000    0.0050   28.7104    5.5000
%10.0000    5.0381    5.0000   75.6615   75.0000    0.5028    0.5000    0.0100   31.0444    8.5000



%%%
figure(1)
hold on
plot(mc_step_ii(1:1000:length(mc_step_ii)), mc_pot_val(1:1000:length(mc_pot_val)), 'b', 'LineWidth', 1.5)
hold off
box on
ylabel('\mbox{Potential value}','Interpreter','latex') % , 'Rotation',0
xlabel('\mbox{MC step}','Interpreter','latex')
%axis([0 n_moves 0 2])
set(gca,'FontSize',16)



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

