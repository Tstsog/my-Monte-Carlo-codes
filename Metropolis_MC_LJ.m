% This Matlab code performs the Monte-Carlo (MC) simulation and computes 
% an average potential energy per particle interacting with the Lennard-Jones
% potential using a Metropolis algorithm. Results for potential energy per particle from this calculation
% are compared with those obtained by Verlet, Ref. [2]. 
%
% Ref. [1] D. Frenkel and B. Smit, "Understanding Molecular Simulation", Acedmic Press (2002);
% Ref. [2] L. Verlet, Phys. Rev. v159, p98 (1967); 
%      
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% July 18, 2024 & University of North Dakota
%
function [] = Metropolis_MC_LJ
clc; close all; 
format short
%
dim_par = 5;        % dimensional parameter 
npart = dim_par^3;  % number of particle 
%
rho_red = 0.880;  % reduced density 
T_red = 0.940;    %  reduced temperature
%
side = (npart/rho_red)^(1/3);  % length of box side
%
sideh = side * 0.5;        
rc = 2.5;                  % cut-off parameter in the Lennard-Jones potential  
rc2 = rc * rc;
%
%dx = 0.100;
dx = side/npart; % an interval, in which particle moves & it can be a small constant (dx = 0.10) as well; & controls an acceptance ratio. 
%
n_iter = 1500.; % number of MC calculation 
%
[x, y, z] = coordinates(dim_par, npart, rho_red);
%
figure(1)
scatter3(x, y, z, "blue", "filled") % initial configuration 
ylabel('y')
xlabel('x')
zlabel('z')
set(gca,'FontSize',16)

% initial potential energy
[epot] = init_pot(x, y, z, side, sideh, rc2, npart);

% Metropolis MC calculation
[pot_value_ave, ratio_accept, x, y, z] = metropolis_lj(x, y, z, npart, n_iter, side, sideh, dx, epot, rc2, T_red);


figure(3)
scatter3(x, y, z, "blue", "filled") % final configuration 
ylabel('y')
xlabel('x')
zlabel('z')
set(gca,'FontSize',16)

%
Potential_energy_tail_correction_per_particle = ((8/3)*pi*rho_red)*((1/3)*(1./rc^9) - (1./rc^3)); % from Ref. [1]. 

pot_value_ave_with_correction = pot_value_ave + Potential_energy_tail_correction_per_particle;

[npart, rho_red, T_red, pot_value_ave, ratio_accept, dx, pot_value_ave_with_correction]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculation potential energy (U*/epsilon) per particle

%%% [npart, rho_red, T_red, pot_value_ave, ratio_accept, dx, pot_value_ave_with_correction]
%
% 125.0000    0.8800    1.0950   -5.2977   73.0411    0.0417   -5.7689 vs -5.66 from Ref. [2]
% 125.0000    0.8800    0.9400   -5.4704   72.3872    0.0417   -5.9416 vs -5.84
% 125.0000    0.8800    0.5910   -5.8751   68.0272    0.0417   -6.3463 vs -6.53
% 125.0000    0.8500    2.2020   -4.3651   79.6235    0.0422   -4.8202 vs -4.76
% 125.0000    0.7500    0.8810   -4.8916   77.0240    0.0440   -5.2932 vs -5.31
% 125.0000    0.6500    0.9000   -4.2394   80.3221    0.0462   -4.5875 vs -4.61
% 125.0000    0.5000    1.3600   -3.0827   86.0619    0.0504   -3.3504 vs -3.38
% 125.0000    0.4500    1.5520   -2.6911   87.3787    0.0522   -2.9320 vs -2.98
% 125.0000    0.3500    1.4180   -2.1014   89.1371    0.0568   -2.2888 vs -2.21 from Ref. [2]

%%%
return
end


%%%
function [x, y, z] = coordinates(dim_par, npart, rho_red)
%
x = zeros(npart,1);
y = zeros(npart,1);
z = zeros(npart,1);
%
side3 = 1./rho_red^(1/3);
%
for i = 1:npart
    j = ceil(i/dim_par^2);
    j1 = j-1;
    k1 = i-j1*dim_par^2;
    k = ceil(k1/dim_par);
    m = rem(i,dim_par) + 1;
      
    x(i,1) = j;
    y(i,1) = k;
    z(i,1) = m;
end
x = x - (dim_par+1)/2;
x = x * side3;
%
y = y - (dim_par+1)/2;
y = y * side3;
%
z = z - (dim_par+1)/2;
z = z * side3;

%%%
return
end

%%%
%%%
function [epot] = init_pot(x, y, z, side, sideh, rc2, npart)
%
epot = 0.;   % potential energy
for i = 1:npart
        %
    for j = i+1:npart
        xx = x(i) - x(j);
        yy = y(i) - y(j);
        zz = z(i) - z(j);
            %
        if (xx <-sideh); xx = xx + side; end %minimum image convention
        if (xx > sideh); xx = xx - side; end 
            %
        if (yy <-sideh); yy = yy + side; end
        if (yy > sideh); yy = yy - side; end             
            %
        if (zz <-sideh); zz = zz + side; end
        if (zz > sideh); zz = zz - side; end 
        %
        r2 = xx * xx + yy * yy + zz * zz;
            %
        if (r2 < rc2)
            %
            r2i = 1/r2;
            r6i = r2i * r2i * r2i;
            %
            epot = epot + 4.*r6i * (r6i - 1.); 
        end %       
    end
end
%
epot = epot/npart; % potential per particle 
%%%
return
end


%%%
function [pot_value_ave, ratio_accept, x, y, z] = metropolis_lj(x, y, z, npart, n_iter, side, sideh, dx, epot, rc2, T_red)
%
sm = 0.;
%
dy = dx;
dz = dx;
%
n_accept= 0.; 
%%%%%%%%%%%%%%%%%%%%%%%%
fileID_save_data_1 = fopen('Metropolis_MC_LJ.txt','w');
%
for ii = 1:n_iter
    %
    for i = 1:npart   %
        
        % old configuration
        x_old = x(i);
        y_old = y(i);
        z_old = z(i);
        
        % potential with old configuration 
        pot_old = 0;

        for j = 1:npart   
           
            if (j ~=i)
                xx = x_old - x(j);
                yy = y_old - y(j);
                zz = z_old - z(j);

                if (xx <-sideh); xx = xx + side; end %minimum image convention
                if (xx > sideh); xx = xx - side; end 
                            %
                if (yy <-sideh); yy = yy + side; end
                if (yy > sideh); yy = yy - side; end             
                            %
                if (zz <-sideh); zz = zz + side; end
                if (zz > sideh); zz = zz - side; end 
                
                r2 = xx * xx + yy * yy + zz * zz;
                 
                if (r2 < rc2)
                    r2i = 1/r2;
                    r6i = r2i * r2i * r2i;
                    pot_old = pot_old + 4 * r6i * (r6i - 1.);
                end
            end
        end   
        %
        x_new = x_old + 2 * dx *(rand(1) - 0.5);
        y_new = y_old + 2 * dy *(rand(1) - 0.5);
        z_new = z_old + 2 * dz *(rand(1) - 0.5);        
      
        %  potential with new configuration 
        pot_new = 0;
        for j = 1:npart
            if (j ~=i)
               
                xx = x_new - x(j);
                yy = y_new - y(j);
                zz = z_new - z(j);

                if (xx <-sideh); xx = xx + side; end %minimum image convention
                if (xx > sideh); xx = xx - side; end 
                            %
                if (yy <-sideh); yy = yy + side; end
                if (yy > sideh); yy = yy - side; end             
                            %
                if (zz <-sideh); zz = zz + side; end
                if (zz > sideh); zz = zz - side; end 
                
                r2 = xx * xx + yy * yy + zz * zz;
    
                if (r2 < rc2)
                    r2i = 1/r2;
                    r6i = r2i * r2i * r2i;
                    pot_new = pot_new + 4 * r6i * (r6i - 1.);
                end
            end
        end
                
        delta_pot = pot_new - pot_old;
        
        if( (delta_pot <= 0) || ( exp(-delta_pot/T_red) > rand(1)) ) % Metropolis scheme 
            %
            sm = sm + delta_pot;
            %
            x(i) = x_new;
            y(i) = y_new;
            z(i) = z_new;
            %
            n_accept = n_accept + 1; 
        end
        %end


    end
    %
    output = [ii,  epot + sm/(npart)];    
%
    fprintf(fileID_save_data_1, '%4.6f \t %8.6f\n', output); 
%
    
end 
fclose(fileID_save_data_1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_mc_data = fopen('Metropolis_MC_LJ.txt', 'r');               % 
read_mc_data = textscan(read_mc_data, '%f %f');
mc_step_ii = read_mc_data{1};
mc_pot_val = read_mc_data{2};

figure(2)
plot(mc_step_ii, mc_pot_val, 'b')
xlabel('\mbox{MC step}','Interpreter','latex') % ,'fontsize',16
ylabel('$U^{\ast}$','Interpreter','latex','Rotation',1) % , 'Rotation',0
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
box on
%
pot_value_ave = sum(mc_pot_val)/length(mc_pot_val);
ratio_accept = 100 * n_accept/(npart*n_iter);  % in percent

%%%
return
end