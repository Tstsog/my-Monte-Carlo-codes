% This Matlab code performs the Monte-Carlo (MC) simulation and computes 
% an average potential energy for a partcile in a double-well potential (V(x)), and
% probability distribution for it in this potential using a Metropolis algorithm [1]. 
%
% Ref. [1] E. Curotto, "Stochastic Simulations of Clusters: Quantum Methods in Flat and Curved Spaces", CPC Press (2010).
%      
%  Double-well potential: V(x) = gamma * x^4 - x^4. 
%
% Average value of potential V(x), V_ave: 
%
%          /
%         | dx*V(x)*exp(-V(x)/(k_B*T))
%         /                                  1
% V_ave = ------------------------------- ~ --- sum_{i}^{N} V(x_{i})
%          /                                 N
%         | dx*exp(-V(x)/(k_B*T))
%         /
%
%
% Probability distribution, f(x):
%
%              exp(-V(x)/(k_B*T))
% f(x) = -------------------------------
%          /
%         | dx*exp(-V(x)/(k_B*T))
%         /
%
%
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% December 30, 2024 & University of North Dakota
%
%
function [] = Metropolis_MC_double_well
%
clc; clear;
format long 
%
gamma = 0.015; % constant in potential V(x)
temp = 3.;     % temperature (i.e. (temp = k_B*T))
nmoves = 10^6; % number of the metropolis moves 
%
% Plot of double well potential V(x)
plot_double_well(gamma);
%
% Obtain an exact probability distribution by computing direct integral form
exact_distribution_double_well_pot(gamma, temp);
%
% Obtain a probability distribution by the metropolis scheme
[Pot_ave, Rejection_ratio] = metropolis_double_well(gamma, temp, nmoves);
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% direct calculation of average of double-well potential V(x)
dx = 0.01;
x = -10:dx:10;
pot = gamma.*x.^4 - x.^2;
sum1 = sum(pot.*exp(-pot./temp).*dx);
sum2 = sum(exp(-pot./temp).*dx);
direct_numerical_calculation_of_pot_average = sum1/sum2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Pot_ave, Rejection_ratio, direct_numerical_calculation_of_pot_average]
%
% [Pot_ave, Rejection_ratio, direct_numerical_calculation_of_pot_average]
% -14.984500894269688  49.160499999999999 -14.992311000773549

%%%
return
end

function [Pot_ave, Rejection_ratio] = metropolis_double_well(gamma, temp, nmoves)
%
deltax = 5.300;
%
x_init = -15.50;
%
nrej = 0.;
sm_pot = 0.;
%
x = x_init;
%
[pot] = double_well(gamma,x);
%
fileID_save_data_1 = fopen('metropolis_double_well_temp3p00.txt','w');
%
for moves = 1:2*nmoves
    %
    xt = x + deltax * (rand(1) - 0.5);
    %
    [pot_t] = double_well(gamma,xt);
    %
    q = exp(-pot_t/temp)/exp(-pot/temp);
    %
    if (rand(1) <= q)
        x = xt ;
        pot = pot_t;
    else
        nrej = nrej + 1;
    end
    %
    if (moves < nmoves)
        %
%        [x];
        output = [moves, x];
        %
        fprintf(fileID_save_data_1, '%4.4f \t %8.12f\n', output);
        %
        sm_pot = sm_pot + pot;
    end
end
%
fclose(fileID_save_data_1);
%
%[sm_pot/nmoves, (nrej/(2*nmoves))*100]

Pot_ave = sm_pot/nmoves;
Rejection_ratio = (nrej/(2*nmoves))*100; % in percent

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_data = fopen('metropolis_double_well_temp3p00.txt', 'r');               % 
read_data = textscan(read_data, '%f %f ');
number_of_moves = read_data{1};
x_val = read_data{2};
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_data = fopen('exact_distribution_double_well_pot.txt', 'r');               % 
read_data = textscan(read_data, '%f %f ');
x_coord = read_data{1};
pot_val = read_data{2};
%
nbins = 30;
figure(2)
hold on
plot(x_coord, pot_val, 'b', LineWidth=1.5)     % exact distribution 
histogram(x_val, nbins,'Normalization','pdf')  % by metropolis scheme
xlabel('$x$','Interpreter','latex') % ,'fontsize',16
ylabel('Probability distribution') % , 'Rotation',0 ,'Rotation',1
hold off
set(gca,'FontSize',18)
box on

%%%
return
end


function exact_distribution_double_well_pot(gamma, temp)
%
x_init = -15.50;
sm = 0.;
dx = 0.01;
nk = 3000;
%
fileID_save_data_1 = fopen('exact_distribution_double_well_pot.txt','w');
%
x = x_init;
%
for k = 1:nk
    %
    [pot] = double_well(gamma,x);
    %
    sm = sm + dx*exp(-pot/temp);
    x = x + dx;
end
%
x = x_init;
for k = 1:nk
    %
    [pot] = double_well(gamma,x);
    %
    fx_dist = exp(-pot/temp)/sm; 
    %
    output = [x, fx_dist];
    %
    x = x + dx;
    %
    fprintf(fileID_save_data_1, '%4.4f \t %8.12f\n', output);
    %
end
%
fclose(fileID_save_data_1);
%%%
return
end


%%%
function plot_double_well(gamma)
%
x_init = -10.0;
dx = 0.05;
nk = 400;
%
fileID_save_data_1 = fopen('plot_double_well.txt','w');
%
x = x_init;
%
for i = 1:nk
    %
    [pot] = double_well(gamma,x);
    %
    output = [x, pot];
    %
    x = x + dx;
    %
    fprintf(fileID_save_data_1, '%4.4f \t %8.12f\n', output);
    %    
end
%
fclose(fileID_save_data_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_data = fopen('plot_double_well.txt', 'r');               % 
read_data = textscan(read_data, '%f %f');
x_coord = read_data{1};
pot_val = read_data{2};
%
figure(1)
hold on
plot(x_coord, pot_val, 'b', LineWidth=1.5)
xlabel('$x$','Interpreter','latex') % ,'fontsize',16
ylabel('$Double\,\,\,well\,\,\, potential$','Interpreter','latex') % , 'Rotation',0 ,'Rotation',1
axis([-10. 10. -20.0 5.00])
hold off
set(gca,'FontSize',18)
box on

%%%
return
end
%

%%%%
function [pot] = double_well(gamma,x)
% Double-well potential
%
pot = gamma.*x^4 - x*x;
%%%
return
end


