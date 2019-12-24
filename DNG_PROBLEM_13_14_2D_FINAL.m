% [Starter code for Problem 8]

% ------------------------------------------------------------
% Solving the bulk with an inhomogeneous iterative method
% ------------------------------------------------------------

clear; clc; close all;

% Parameters:
z_nn = 4;   % Number of nearest-neighbour in lattice (square = 4).
z_nnn = 4;  % Number of next-nearest-neighbours in lattice (square = 4).
Lx = 102;    % Number of sites along x-axis.
Ly = 42;    % Number of sites along y-axis.
sigma = 1;  % Size of a site (defines our units of length).
beta = 1.2; % Inverse temperature beta*epsilon.
mu = -2.5; % Chemical potential mu/epsilon.
J = linspace(1, 11, 11);%J points for the line graph plot
bew = 1.3;

potential = zeros(Ly);
attract = 1.6; %wall attraction constant
k = 1;         %wall depth

rho_0 = 0.4;   % Initial density.
tol = 1e-10;   % Convergence tolerance.
count = 100000; % Upper limit for iterations.
alpha  = 0.03; % Mixing parameter.
Dx = 20;
Dy = 30;
N = Dx*Dy;

conv = 1; cnt = 1;       % Convergence value and counter.
rho = zeros(Lx,Ly); % Initialise rho to the starting guess(i-th rho_old) in Eq(47)
rngx = (Lx/2 - Dx/2):(Lx/2 + Dx/2 - 1);
rngy = 2:(Dy+1);
rho(rngx,rngy) = 1;
rho_rhs = zeros(Lx,Ly);  % Initialise rho_new to zeros.



% Solve equations iteratively:
while conv>=tol && cnt<count
  cnt = cnt + 1; % Increment counter.
  % Loop over all lattice sites:
  rho_prev = rho;
  for i=1:Lx
    for j=2:(Ly-1)
        %Defining the Lennard-Jones potential
        potential(j) = -bew*(j-1)^(-3); 
       
        % Handle the periodic boundaries for x and y:
        left = mod((i-1)-1,Lx) + 1;  % i-1, maps 0 to Lx.
        right = mod((i+1)-1,Lx) + 1; % i+1, maps Lx+1 to 1.

        rho_rhs(i,j) = (1 - rho(i,j))*exp((beta*(rho(i,j-1) + rho(i,j+1) + rho(left,j) + rho(right,j) + (1/4)*(rho(left,j-1) + rho(right,j-1) + rho(left,j+1) + rho(right,j+1))+mu) - potential(j)));
    end
  end
  
  
    rho = alpha*rho_rhs + (1 - alpha)*rho; % Mix the new and current solutions for next iteration.
    rho = rho*N/(sum(sum(rho)));
    conv = sum(sum((rho - rho_prev).^2)); % Convergence value is the sum of the differences between new and current solution.
   
    if mod(cnt,1000)==0
        disp(['conv = ' num2str(conv) ' cnt = ' num2str(cnt)]); % Display final answer.
    end
end
x=[1:100]
y=2.1445*x - 127.5
figure(1);
pcolor(0:(Lx-1),1:(Ly-1), rho(:,2:end)');
shading interp;
axis equal;
ylim([1,Ly-1]);
xlabel('x/\sigma','FontSize', 18); ylabel('y/\sigma','FontSize', 18);
c = colorbar;
c.Label.String = '\rho\sigma^2';
c.Label.FontSize = 18;
title('Droplet shape', 'FontSize', 15)
% title('Initial droplet shape', 'FontSize', 15)
txt = ['\beta\epsilon_{w} = ' num2str(bew) ''];
text(10,30,txt,'color','w','HorizontalAlignment','left','FontSize', 20)


hold on;
contour(0:(Lx-1),1:(Ly-1),rho(:,2:end)',[0.5,0.5],'r-','LineWidth',1.5)
% plot(x,y,'k', 'LineWidth', 2.5)
hold off;

% run to line 62 then type
% pcolor(rho)
% pcolor(rho')
% shading interp;
% shg