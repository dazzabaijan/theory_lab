% [Starter code for Problems 6 and 7]

% ------------------------------------------------------------
% Plotting bulk homogeneous liquid grand potential minimum
% ------------------------------------------------------------

clear; clc; close all;

% Parameters:
sigma = 1;                % Size of a site (defines our units of length).
z_nn = 4;                 % Number of nearest-neighbour in lattice (square = 4).
z_nnn = 4;                % Number of next-nearest-neighbours in lattice (square = 4).
samp = 100;               % Sampling of the density.
N = 50;
rhox = -1.4279*ones(N);
rhoy = linspace(0,1,N);
spinodalx = -2.5*ones(N);
spinodaly = linspace(0,1,N);
rho = linspace(0,1,samp); % The range of allowed densities.
beta = 1/0.4 % Inverse temperature beta*epsilon.(epsilon/Kb*T)
crit_beta = (5/4)*ones(N) 
mu = linspace(-5,0.001,N); % Chemical potential mu/epsilon.
p1 = 0;
p2 = 0;

%App = approximate rho
app_high = zeros(N);
app_mid = zeros(N);
app_low = zeros(N);
app_diff = zeros(N);
mu_plus = linspace(0,N,N)
mu_minus = linspace(0,N,N)
mu_coex = zeros(N);


for j = 1:N

    % Create a Matlab function handle f for the function bulk_func() which is the
    % lefthand side of Eq. (45) giving the bulk grand potential minimum:
    f = @(x) bulk_func(x,beta,mu(j),z_nn,z_nnn); % Here x is the variable argument.

    % Solve this f = 0 using fzero():

    %Low density gas
    %Solution coming from a low-density initial guess:
    app_low(j) = fzero(f,1e-4); 

    %Solution coming from an equal density initial guess: 
    app_mid(j) = fzero(f,0.5);

    %High density liquid
    %Solution coming from a high-density initial guess:
    app_high(j) = fzero(f,1 - 1e-4);

    %Mix both graphs
    app_diff(j) = app_high(j) - app_low(j);
    if app_diff(j) >0.01 
        app_diff(j) = nan;

    else app_diff(j) = app_low(j);
    end
end





%At each point of every beta values, rho+&rho- points can be calculated
rho_plus = (1/10)*(5 + ((5)^(0.5))*(5-(4/beta))^(0.5))
rho_minus = (1/10)*(5 - ((5)^(0.5))*(5-(4/beta))^(0.5))

%At each point of every rho+&rho- points(and hence beta points),
%mu+&mu- can be calculated
mu_plus = (1/beta)*log(rho_plus/(1-rho_plus)) - 5*rho_plus
mu_minus = (1/beta)*log(rho_minus/(1-rho_minus)) - 5*rho_minus
    

app_high(app_high == 0) = NaN;
app_mid(app_mid == 0) = NaN;
app_low(app_low == 0) = NaN;

% figure (4)
% plot(mu, app_mid, 'LineWidth', 3)
% xlim([-5,0]);
% 
% figure (5)
% plot(mu, app_low, 'b--', 'LineWidth',3)
% xlim([-5,0]);
% 
% figure (6)
% plot(mu, app_mid, 'g--', 'LineWidth',3)
% xlim([-5,0]);

figure (7)
hold on;
plot(mu, app_high,'r-', 'LineWidth',2)
plot(mu, app_high,'o','MarkerSize',5)
plot(mu, app_low, 'b-', 'LineWidth',2)
plot(mu, app_low, 'x','MarkerSize',5)
plot(mu, app_mid, 'g-', 'LineWidth',2)
plot(mu, app_mid, 'd','MarkerSize',5)
% plot(rhox, rhoy, '--', 'LineWidth',2)
plot(spinodalx, spinodaly, '--', 'LineWidth',3);
xlim([-5,0]);
title(['Bulk phase diagram with solutions for \beta\epsilon = ' num2str(1/beta)])
xlabel('\mu(\rho)/\epsilon','FontSize', 18); ylabel('\rho','FontSize', 18);
h = zeros(4, 1);
h(1) = plot(NaN,NaN,'--black','linewidth',2);
h(2) = plot(NaN,NaN,'-r', 'linewidth',2);
h(3) = plot(NaN,NaN,'-g', 'linewidth',2);
h(4) = plot(NaN,NaN,'-b', 'linewidth',2);
legend(h, 'Binodal','Low \rho solution','Equal \rho solution','High \rho solution', 'Location','northwest');
hold off;
txt_2 = {'gaseous','phase'};
text(-4.2,0.5,txt_2,'HorizontalAlignment','center','FontSize', 15, 'color', 'b')

txt_3 = {'liquid','phase'};
text(-0.8,0.5,txt_3,'HorizontalAlignment','center','FontSize', 15, 'color', 'b')

% an = annotation('doublearrow');
% an.X = [-4 -3.5]
% an.Y = [0.05 0.05]
% c = an.Color;
% an.Color = 'red';

figure(8)
hold on;
plot(app_high, mu, 'r-', 'LineWidth',2)
plot(app_high, mu, 'o','MarkerSize',5)
plot(app_low, mu,  'b-', 'LineWidth',2)
plot(app_low, mu, 'x','MarkerSize',5)
plot(app_mid, mu,  'g-', 'LineWidth',2)
plot(app_mid, mu, 'd','MarkerSize',5)
plot(rhox, rhoy, '--', 'LineWidth',2)
hold off;
ylim([-5,0]);
xlim([0,1]);
title(['Bulk phase diagram with high, mid and low density solutions for \beta\epsilon = ' num2str(1/beta)])
ylabel('\mu(\rho)/\epsilon','FontSize', 18); xlabel('\rho','FontSize', 18);


% figure(8)
% plot(mu, app_diff, 'LineWidth',3)
% xlim([-5,0]);

% % Now evaluate the f over density:
% f_vals = f(rho);
% s
% % Plot f(rho) and the three (not necessarily unique) solutions found above:
% figure(1);
% set(gcf,'Color','white')
% plot(rho,f_vals,'r-',rho,zeros(1,samp),'k--',rho_low,0,'bd',rho_mid,0,'go',rho_high,0,'mx','LineWidth',3,'MarkerSize',15)
% set(gca,'FontSize',14)
% legend({['\mu/\epsilon = ' num2str(mu)]},'Location','northwest')
% title(['Grand potential minimum with \rho for \beta\epsilon = ' num2str(beta)])
% xlabel('\rho'); ylabel('f(\rho)');