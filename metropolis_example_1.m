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
n_moves = 500000; % number of Metropolis moves
delta_x = 3.50;   % parameter used to control a rejection rate between 40 and 60%
%
x = 0.0;          % initial coordinate
n_rej = 0.;       % number of rejection 
sm_v = 0.;        % sum of V
sm_v2 = 0.;       % sum of V^2
%
%%%%%%%%%%%%%%%%%%%%%%%%
fileID_save_data_1 = fopen('metropolis_example_1.txt','w');
%
x4_save = zeros(n_moves, 1);
x4_sq_save = zeros(n_moves, 1);
%
for ii = 1:n_moves
    %
    xt = x + delta_x * (rand(1) - 0.5);
    %
    q = exp(-xt^2)/exp(-x^2); % probability ratio, q = exp(-V(xt))/exp(-V(x))
    if (rand(1) < q)
        x = xt;
        %
        sm_v = x^4;
        sm_v2 = x^4 * x^4;        
    else
        n_rej = n_rej + 1;
    end
    %


    %%%
    x4_save(ii) = sm_v;
    x4_sq_save(ii) = sm_v2;
    output = [ii, x4_save(ii), x4_sq_save(ii)];
    %
    fprintf(fileID_save_data_1, '%4.6f \t %4.6f \t %8.6f\n', output); 
end
fclose(fileID_save_data_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_mc_data = fopen('metropolis_example_1.txt', 'r');               % 
read_mc_data = textscan(read_mc_data, '%f %f %f');
mc_step_ii = read_mc_data{1};
mc_x2_val = read_mc_data{2};
mc_x4_val = read_mc_data{3};
%

%%%
x4_ave = sum(mc_x2_val)/length(mc_x2_val);
x4_sq_ave = sum(mc_x4_val)/length(mc_x4_val);
sigma_std = sqrt((x4_sq_ave - x4_ave^2)/length(mc_x2_val));

%%%
x4_exact = 0.75;                                  % exact value  
rejection_ratio = 100*(n_rej/n_moves);            % rejection in percent
%
[x4_ave, x4_exact, sigma_std, rejection_ratio ]
%
%[x4_ave, x4_exact, sigma_std, rejection_ratio ]
% 0.7497    0.7500    0.0035   43.8568

%%%%%%%%%%%%%

figure(1)
hold on
plot(mc_step_ii(1:1000:length(mc_step_ii)), mc_x2_val(1:1000:length(mc_x2_val)), 'b')
hold off
box on
ylabel('$x^{2}$ \mbox{value}','Interpreter','latex') % , 'Rotation',0
xlabel('\mbox{MC step}','Interpreter','latex')
%axis([0 n_moves 0 2])
set(gca,'FontSize',16)



%%%
return
end



