% [Starter code for Problem 8]

% ------------------------------------------------------------
% Solving the bulk with an inhomogeneous iterative method
% ------------------------------------------------------------

clear; clc; close all;

% Parameters:
z_nn = 4;   % Number of nearest-neighbour in lattice (square = 4).
z_nnn = 4;  % Number of next-nearest-neighbours in lattice (square = 4).
Lx = 4;    % Number of sites along x-axis. Was 40
Ly = 20;    % Number of sites along y-axis. Was 40
sigma = 1;  % Size of a site (defines our units of length).
mu_coex = -2.5;
beta = 1.2; % Inverse temperature beta*epsilon.
% mu = [-2.67, -2.6, -2.53, -2.5];     % Chemical potential mu/epsilon.
marker = ['+', '*', 'o', 'd', 'x', '^']
% mu = mu_coex - [0.17, 0.104, 0.03, 0].*beta;
mu = mu_coex - [0.2, 0.17, 0.10, 0.09, 0.03, 0];
% mu = mu_coex - [0.2,0.104,0.1,0.04,0.004,0].*beta; % Chemical potential mu/epsilon.
% marker = ['+', '*', 'o', 'd', 'x', '^']
% mu = [-2.53]
% marker = ['+']
J = linspace(1, 10, 10);             %J points for the line graph plot

potential = zeros(Ly);
gamma = zeros(Ly, length(mu));
bew = 1.6;                      %wall attraction constant
k = 1;                              %wall depth

rho_0 = 0.4;                        % Initial density.
tol = 1e-12;                        % Convergence tolerance.
count = 30000;                      % Upper limit for iterations.
alpha = 0.01;                       % Mixing parameter.

rho = rho_0*ones(Ly); % Initialise rho to the starting guess(i-th rho_old) in Eq(47)
rho_rhs = zeros(Ly);  % Initialise rho_new to zeros.

figure(1);
hold on;
for i=1:length(mu)
    conv = 1; cnt = 1; % Convergence value and counter.
    
    %Solve equations iteratively:
    while conv>=tol && cnt<count
        cnt = cnt + 1; % Increment counter.
        %Loop over all lattice sites:
        for j=1:Ly
            %Defining the Lennard-Jones potential
            if j<k
                potential(j) = 1000000000;
            else
                potential(j) = -bew*(j-k).^(-3); 
            end
            %Handle the periodic boundaries for x and y:
            %left = mod((i-1)-1,Lx) + 1; % i-1, maps 0 to Lx.
            %right = mod((i+1)-1,Lx) + 1; % i+1, maps Lx+1 to 1.
            if j<k+1 %depth of wall
                rho_rhs(j) = 0;
                rho(j) = 0;
            elseif j<(20+k)
                rho_rhs(j) = (1 - rho(j))*exp((beta*((3/2)*rho(j-1) + (3/2)*rho(j+1) + 2*rho(j) + mu(i)) - potential(j)));
            else
                rho_rhs(j) = rho_rhs(j-1);
            end
%             disp(j);
        end
            
        conv = sum(sum((rho - rho_rhs).^2));   % Convergence value is the sum of the differences between new and current solution.
        rho = alpha*rho_rhs + (1 - alpha)*rho; % Mix the new and current solutions for next iteration.
       
    end
%     disp(i);
    
    disp(['conv = ' num2str(conv) ' cnt = ' num2str(cnt)]); % Display final answer.
    plot(rho(1:10), 'LineStyle', '-', 'Marker', marker(i), 'LineWidth', 1.5);
end
hold off;
xlabel('y/\sigma','FontSize', 18); ylabel('\rho\sigma^2','FontSize', 18);

leg = legend('-2.6667', '-2.6000', '-2.5333', '-2.5000' );
leg = legend('-0.108', '-0.104', '-0.098', '-0.094', '-0.03', '0' );
% legend({['\mu/\epsilon = ' num2str(mu) ,' , ', '\beta\epsilon = ' num2str(beta)]});
% legend({['\mu/\epsilon = ' num2str(mu) ,' , ', '\beta\epsilon = ' num2str(beta)]},'Location','northwest');

% leg = legend('-0.2', '-0.104', '-0.1', '-0.04', '0.004', '0');
% title(leg, {['\beta\epsilon = ' num2str(beta)], ['\beta*(\mu - \mu_{coex}) = ']});
% title(leg, {['\beta\epsilon = ' num2str(beta)], ['\mu/\epsilon = ']});
title(leg, {['\beta\epsilon = ' num2str(beta)], ['\beta*(\mu - \mu_{coex})']});

% title(leg, {['\beta\epsilon = ' num2str(beta)]});
% title('Density profiles close to the wall(1D)','FontSize', 15)
title({['Density profiles close to the wall with adsorption \beta\epsilon_w = ' num2str(bew)]},'FontSize', 12)
% xlim([1 inf])

% figure(2)
% pcolor(rho');
% % title('Density profiles close to the wall(1D)','FontSize', 15)
% xlabel('y/\sigma','FontSize', 18); ylabel('x/\sigma','FontSize', 18);
% h = colorbar;
% h.Label.String = '\rho\sigma^2';
% h.Label.FontSize = 18;
% ylim([1 4])
% title({['Density profiles close to the wall(2D) \beta\epsilon_w = ' num2str(bew)]},'FontSize', 12)
