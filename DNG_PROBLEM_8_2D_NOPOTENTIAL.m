% [Starter code for Problem 8]

% ------------------------------------------------------------
% Solving the bulk with an inhomogeneous iterative method
% ------------------------------------------------------------

clear; clc; close all;

% Parameters:
z_nn = 4;   % Number of nearest-neighbour in lattice (square = 4).
z_nnn = 4;  % Number of next-nearest-neighbours in lattice (square = 4).
Lx = 40;    % Number of sites along x-axis.
Ly = 40;    % Number of sites along y-axis.
sigma = 1;  % Size of a site (defines our units of length).
beta = 1; % Inverse temperature beta*epsilon.
mu = [-2]; % Chemical potential mu/epsilon.
% mu = [-2.67, -2.6, -2.53, -2.5]
marker = ['+']
% marker = ['+', '*', 'o', 'd']
J = linspace(1, 11, 11);%J points for the line graph plot

attract = 1.6; %wall attraction constant
k = 1;         %wall depth


rho_0 = 0.4;   % Initial density.
tol = 1e-12;   % Convergence tolerance.
count = 30000; % Upper limit for iterations.
alpha  = 0.01; % Mixing parameter.


rho = rho_0*ones(Lx,Ly); % Initialise rho to the starting guess(i-th rho_old) in Eq(47)
rho_rhs = zeros(Lx,Ly);  % Initialise rho_new to zeros.


figure(1);
subplot(1,2,1);
hold on;

for loop = 1:length(mu)
    conv = 1; cnt = 1;       % Convergence value and counter.
    % Solve equations iteratively:
    while conv>=tol && cnt<count
        cnt = cnt + 1; % Increment counter.
        % Loop over all lattice sites:
        for i=1:Lx
            for j=1:Ly
                %Defining the Lennard-Jones potential
                
                % Handle the periodic boundaries for x and y:
                left = mod((i-1)-1,Lx) + 1; % i-1, maps 0 to Lx.
                right = mod((i+1)-1,Lx) + 1; % i+1, maps Lx+1 to 1.
                if j<k+1 %depth of wall
                    rho_rhs(i,j) = 0;
                    rho(i,j) = 0;
                elseif j<(20+k)
                    rho_rhs(i,j) = (1 - rho(i,j))*exp((beta*(rho(i,j-1) + rho(i,j+1) + rho(left,j) + rho(right,j) + (1/4)*(rho(left,j-1) + rho(right,j-1) + rho(left,j+1) + rho(right,j+1)) + mu(loop))));
                else
                    rho_rhs(i,j) = rho_rhs(i,j-1);
                end
            end
        end

        conv = sum(sum((rho - rho_rhs).^2)); % Convergence value is the sum of the differences between new and current solution.
        rho = alpha*rho_rhs + (1 - alpha)*rho; % Mix the new and current solutions for next iteration.
    end
    disp(['conv = ' num2str(conv) ' cnt = ' num2str(cnt)]); % Display final answer.
    plot(J, rho(1,1:11), 'LineStyle', '-', 'Marker', marker(loop), 'LineWidth', 1.5);
end
hold off;
legend({['\mu/\epsilon = ' num2str(mu) ,' , ', '\beta\epsilon = ' num2str(beta)]},'Location','northwest');
xlabel('y/\sigma','FontSize', 18); ylabel('\rho\sigma^2','FontSize', 18);
title('Density profiles close to the wall(2D model)','FontSize', 15)
xlim([1 inf])

subplot(1,2,2);
pcolor(rho);
title('Density profiles close to the wall(2D model)','FontSize', 15)
xlabel('y/\sigma','FontSize', 18); ylabel('x/\sigma','FontSize', 18);
h = colorbar;
h.Label.String = '\rho\sigma^2';
h.Label.FontSize = 20;

% figure(2);
% pcolor(rho);
% title('Density profiles close to the wall(2D model)','FontSize', 15)
% xlabel('y/\sigma','FontSize', 18); ylabel('x/\sigma','FontSize', 18);
% h = colorbar;
% h.Label.String = '\rho\sigma^2';
% h.Label.FontSize = 20;
% 

% For a different lattice size, change Lx = n, Ly = m, line 20 size of array J, line 68 rho(1,1:n), line 56 j<((m-1)+k)
