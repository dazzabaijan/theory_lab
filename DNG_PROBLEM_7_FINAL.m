% [Starter code for Problems 6 and 7]

% ------------------------------------------------------------
% Plotting bulk homogeneous liquid grand potential minimum
% ------------------------------------------------------------

clear; clc; close all;

% Parameters:
sigma = 1;                % Size of a site (defines our units of length).
z_nn = 4;                 % Number of nearest-neighbour in lattice (square = 4).
z_nnn = 4;                % Number of next-nearest-neighbours in lattice (square = 4).
samp = 500;               % Sampling of the density.
N = 50;
rho = linspace(0,1,samp); % The range of allowed densities.
beta = 1./linspace(1e-2,2,N) % Inverse temperature beta*epsilon.(epsilon/Kb*T)
crit_beta = (5/4)*ones(N) 
mu = linspace(-5,0.001,N); % Chemical potential mu/epsilon.
p1 = 0;
p2 = 0;

%App = approximate rho
app_high = zeros(N,N);
app_mid = zeros(N,N);
app_low = zeros(N,N);
app_diff = zeros(N,N);
mu_plus = linspace(0,N,N)
mu_minus = linspace(0,N,N)
mu_coex = zeros(N);


for i = 1:N
    for j = 1:N
    
        % Create a Matlab function handle f for the function bulk_func() which is the
        % lefthand side of Eq. (45) giving the bulk grand potential minimum:
        f = @(x) bulk_func(x,beta(i),mu(j),z_nn,z_nnn); % Here x is the variable argument.

        % Solve this f = 0 using fzero():
           
        %Low density gas
        %Solution coming from a low-density initial guess:
        app_low(i,j) = fzero(f,1e-4); 

        %Solution coming from an equal density initial guess: 
        app_mid(i,j) = fzero(f,0.5);
        
        %High density liquid
        %Solution coming from a high-density initial guess:
        app_high(i,j) = fzero(f,1 - 1e-4);
        
        %Mix both graphs
        app_diff(i,j) = app_high(i,j) - app_low(i,j);
        if app_diff(i,j) >0.01 
            app_diff(i,j) = nan;
            
        else app_diff(i,j) = app_low(i,j);
        end
    end
end



for i = 1:N
    %At each point of every beta values, rho+&rho- points can be calculated
    rho_plus = (1/10)*(5 + ((5)^(0.5))*(5-(4/beta(i)))^(0.5))
    rho_minus = (1/10)*(5 - ((5)^(0.5))*(5-(4/beta(i)))^(0.5))
    
    %At each point of every rho+&rho- points(and hence beta points),
    %mu+&mu- can be calculated
    mu_plus(i) = (1/beta(i))*log(rho_plus/(1-rho_plus)) - 5*rho_plus
    mu_minus(i) = (1/beta(i))*log(rho_minus/(1-rho_minus)) - 5*rho_minus
    
    %mu+ and mu- solutions are imaginary at the apex
    %if the absolute value of mu+&mu- have an imaginary part
    %(i.e. they ARE imaginary, then they're set to none)
    if abs(imag(mu_plus(i))) > 0
        mu_plus(i) = nan;
        
        if p1==0
            %find out the i-th term of the mu+ for which it's imaginary
            %then set that as root1
            root1 = i;
            %set pl=0 now that the lower bound of i-th value for which mu+
            %is 0 is found so that it won't get replaced again if the loop
            %keeps running up to N-loops
            p1 = 1;
        end
    end
    
    if abs(imag(mu_minus(i))) > 0
        mu_minus(i) = nan;
        if p2==0
            %find out the i-th term of the mu- for which it's imaginary
            %then set that as root2
            root2 = i;
            %set pl=0 now that the lower bound of i-th value for which mu-
            %is 0 is found so that it won't get replaced again if the loop
            %keeps running up to N-loops
            p2 = 1;
        end
    end
   
end

for i=1:N
    if i > root1
        mu_coex(i) = nan;

    elseif i > root2
        mu_coex(i) = nan;

    else 
        mu_coex(i) = -5/2;

    end

end

%set the elements of mu_coex that's 0 into NaN otherwise there's a weird
%line plotted at mu_coex = 0 (somehow it is also a solution)
mu_coex(mu_coex == 0) = NaN;

% figure(1)
% subplot(2,1,1)
% pcolor(mu, 1./beta, app_high); shading interp;
% c = colorbar;
% c.Label.String = '\rho\sigma^2';
% c.Label.FontSize = 20;
% xlabel('\mu/\epsilon','FontSize', 18); ylabel('k_{B}T/\epsilon','FontSize', 18);
% title('The bulk phase diagram(high density initial guess)','FontSize', 15)
% txt = {'\rho_0 = 0.9999'};
% text(-0.3,1.7,txt,'HorizontalAlignment','right','FontSize', 11)
% 
% subplot(2,1,2)
% pcolor(mu, 1./beta, app_low); shading interp;
% c = colorbar;
% c.Label.String = '\rho\sigma^2';
% c.Label.FontSize = 20;
% xlabel('\mu/\epsilon','FontSize', 18); ylabel('k_{B}T/\epsilon','FontSize', 18);
% title('The bulk phase diagram(low density initial guess)','FontSize', 15)
% txt = {'\rho_0 = 0.0001'};
% text(-0.3,1.7,txt,'HorizontalAlignment','right','FontSize', 11)
% 
% figure(2)
% pcolor(mu, 1./beta, app_mid); shading interp;
% c = colorbar;
% c.Label.String = '\rho\sigma^2';
% c.Label.FontSize = 20;
% xlabel('\mu/\epsilon','FontSize', 18); ylabel('k_{B}T/\epsilon','FontSize', 18);
% title('The bulk phase diagram(equal density initial guess)','FontSize', 15)
% txt = {'\rho_0 = 0.5'};
% text(-0.4,1.7,txt,'HorizontalAlignment','right','FontSize', 11)
% 
% figure(3)
% pcolor(mu,1./beta, app_diff); shading interp;
% c = colorbar;
% c.Label.String = '\rho\sigma^2';
% c.Label.FontSize = 20;
% 
% hold on;
% plot(mu_plus, 1./beta, 'Linewidth', 4, 'color', 'r')     %Plot the evaluated mu+ values for each loop
% plot(mu_minus, 1./beta, 'Linewidth', 4, 'color', 'r')    %Plot the evaluated mu- values for each loop
% plot(mu_coex, 1./beta, 'color', 'black', 'linewidth', 2, 'linestyle', '--');
% 
% %Plot the horizontal critical temperature line
% plot(mu_plus, crit_beta, 'color', 'm', 'linewidth', 2, 'linestyle', '--');
% plot(mu_minus, crit_beta, 'color', 'm', 'linewidth', 2, 'linestyle', '--');
% 
% h = zeros(3, 1);
% h(1) = plot(NaN,NaN,'-r','linewidth',2);
% h(2) = plot(NaN,NaN,'--black', 'linewidth',2);
% h(3) = plot(NaN,NaN,'--m', 'linewidth',2);
% legend(h, 'Spinodal','Binodal \mu/\epsilon = -2.5','T_{c} = 1.25');
% hold off;
% 
% txt = {'\mu_{coex} \rightarrow'};
% text(-5/2,0.3,txt,'HorizontalAlignment','right','FontSize', 15)
% 
% txt_2 = {'gaseous','phase'};
% text(-4,0.8,txt_2,'HorizontalAlignment','center','FontSize', 15, 'color', 'w')
% 
% txt_3 = {'liquid','phase'};
% text(-1,0.8,txt_3,'HorizontalAlignment','center','FontSize', 15, 'color', 'b')
% 
% txt_4 = {'unstable','phase-','separated','region'};
% text(-5/2,0.5,txt_4,'HorizontalAlignment','center','FontSize', 10, 'color', 'black')
% 
% txt_5 = {'Critical temperature', '\downarrow'};
% text(-1,1.35,txt_5,'HorizontalAlignment','center','FontSize', 11)
% 
% title('The bulk phase diagram for the 2D lattice fluid','FontSize', 15)
% xlabel('\mu/\epsilon','FontSize', 18); ylabel('k_{B}T/\epsilon','FontSize', 18);

figure(4);
plot(rho, 1./beta);




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