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
beta = 1.2; % Inverse temperature beta*epsilon.
mu = -2.67; % Chemical potential mu/epsilon.
J = linspace(1, 11, 11);%J points for the line graph plot

potential = zeros(Ly);
attract = 1.6; %wall attraction constant
k = 1;         %wall depth


rho_0 = 0.4;   % Initial density.
tol = 1e-12;   % Convergence tolerance.
count = 30000; % Upper limit for iterations.
alpha  = 0.01; % Mixing parameter.
 
conv = 1; cnt = 1;       % Convergence value and counter.
rho = rho_0*ones(Lx,Ly); % Initialise rho to the starting guess(i-th rho_old) in Eq(47)
rho_rhs = zeros(Lx,Ly);  % Initialise rho_new to zeros.

% Solve equations iteratively:
while conv>=tol && cnt<count
  cnt = cnt + 1; % Increment counter.
  % Loop over all lattice sites:
  for i=1:Lx
    for j=1:Ly
        %Defining the Lennard-Jones potential
        if j<k
            potential(j) = 1000000000
        else
            potential(j) = -attract*(j-k)^(-3); 
        end
      % Handle the periodic boundaries for x and y:
      left = mod((i-1)-1,Lx) + 1; % i-1, maps 0 to Lx.
      right = mod((i+1)-1,Lx) + 1; % i+1, maps Lx+1 to 1.
      if j<k+1 %depth of wall
          rho_rhs(i,j) = 0;
          rho(i,j) = 0;
      elseif j<(20+k)
          rho_rhs(i,j) = (1 - rho(i,j))*exp((beta*(rho(i,j-1) + rho(i,j+1) + rho(left,j) + rho(right,j) + (1/4)*(rho(left,j-1) + rho(right,j-1) + rho(left,j+1) + rho(right,j+1)) + mu) - potential(j)));
      else
          rho_rhs(i,j) = rho_rhs(i,j-1);
      end
    end
  end
  
  conv = sum(sum((rho - rho_rhs).^2)); % Convergence value is the sum of the differences between new and current solution.
  rho = alpha*rho_rhs + (1 - alpha)*rho; % Mix the new and current solutions for next iteration.
end



disp(['conv = ' num2str(conv) ' cnt = ' num2str(cnt)]); % Display final answer.
figure(1);
pcolor(rho);

figure(2);
plot(J, rho(1,1:11));
hold on;
mu = -2.53
plot(J, rho(1,1:11));
hold off;