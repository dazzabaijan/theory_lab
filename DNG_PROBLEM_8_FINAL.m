% [Starter code for Problem 8]

% ------------------------------------------------------------
% Solving the bulk with an inhomogeneous iterative method
% ------------------------------------------------------------

clear; clc; close all;

% Parameters:
z_nn = 4;   % Number of nearest-neighbour in lattice (square = 4).
z_nnn = 4;  % Number of next-nearest-neighbours in lattice (square = 4).
Lx = 4;    % Number of sites along x-axis.
Ly = 20;    % Number of sites along y-axis.
sigma = 1;  % Size of a site (defines our units of length).
beta = 1.2; % Inverse temperature beta*epsilon.
mu = [-2.67, -2.6, -2.53, -2.5];     % Chemical potential mu/epsilon.
marker = ['+', '*', 'o', 'd']
J = linspace(1, 10, 10);             %J points for the line graph plot

potential = zeros(Ly);
gamma = zeros(Ly, length(mu));
attract = 1.6;                      %wall attraction constant
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
                potential(j) = -attract*(j-k).^(-3); 
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
legend('\mu/\epsilon = -2.6667', '\mu/\epsilon = -2.6000', '\mu/\epsilon = -2.5333', '\mu/\epsilon = -2.5000' );
xlabel('y/\sigma','FontSize', 18); ylabel('\rho\sigma^2','FontSize', 18);
title('Density profiles close to the wall','FontSize', 15)
figure(2)
pcolor(rho);


