close all
clear 
clc

%==========================================================================
% MOCAT 3: propagation 

% To compare with workbooks:
% - set same data/species properties
% - set initial conditions on S,D,N,lambda

% To check:
% Ask Miles below for d_upper/n_upper at k=N_shell

%==========================================================================
% Data
%==========================================================================

N_shell = 5; % with gradient use: 20, 70, 140
h_min = 200; % initial altitude [km]
h_max = 900; % final altitude [km] 2000
R0 = linspace(h_min,h_max,N_shell+1);
R02 = linspace(h_min,h_max,N_shell+1);
deltaH = R02(2)-R02(1); % thickness of the shell [km]
re = 6378.1; % [km]
years_s = 365*24*3600;

mu = 3.986004418e14;%meters^3/s^2
R1 = R0;
R0 = (re+R0)*1000; % [m]
V = 4/3*pi*(diff(R0.^3)); % volume of the shells [m^3]

Dt = 5; % operation lifetime: after those years, active satellites become derelict [years] 8
v_imp = 10; % impact velocity [km/s]
v = v_imp*1000*(24*3600*365.25);% impact velocity [m/year]
alpha = 0.2; % fraction of collisions that satellites fail to avoid 0.2
alpha_active = 0.01; % fraction of collisions occuring between active satellites 0.01
delta = 10; %10
P = 0.95;
Dhl = deltaH*1000;
Dhu = -deltaH*1000;
Cd = 2.2;

%==========================================================================
% physical characteristics
%==========================================================================

mass      = [223 223 0.640];      % median mass values from MASTER [kg] SOMMA
A         = [1.741 1.741 0.020]; % median area values from MASTER [m^2] SOMMA
diameter  = [1.490 1.490 0.180]; % median diameter values from MASTER [m] SOMMA

n_obj = length(mass);
radii = diameter/2;
area_mass = A./mass; % area/mass ratio  [# shells, # object types]
beta = Cd*area_mass; % ballistic coefficient [# shells, # object types]

%==========================================================================
% NASA standard breakup models implementation for collisions
%==========================================================================

LC = 0.1; % minimum size of fragments [m]
n_f_catastrophic = @(M1,M2) 0.1*LC^(-1.71)*(M1+M2)^(0.75); % number of fragments generated during a catastrophic collision (NASA standard break-up model). M is the sum of the mass of the objects colliding in kg
n_f_damaging = @(M1,M2) 0.1*LC^(-1.71)*(min(M1,M2)*v_imp^2)^(0.75); % number of fragments generated during a non-catastrophic collision (improved NASA standard break-up model: takes into account the kinetic energy). M is the mass of the less massive object colliding in kg
K0 = zeros(n_obj,n_obj); % number of fragments generated in each collisions among the species for each altitude bin
ind_obj_catastrophic = 2; % hypothesis of catastrophic collisions only between satellites and derelict satellites
sigma = zeros(n_obj,n_obj);
for q = 1:length(radii)
    mass_q = mass(q);
    for qq = 1:length(radii)
        sigma(q,qq) = (radii(q)+radii(qq))^2;
        mass_qq = mass(qq);
        if q<=ind_obj_catastrophic && qq<=ind_obj_catastrophic
            K0(q,qq) = n_f_catastrophic(mass_q,mass_qq); % catastrophic collision
        else
            K0(q,qq) = n_f_damaging(mass_q,mass_qq); % (non-catastrophic) damaging collision
        end
    end
end

% phi=sigma*v; % adopted before
f_corrective = 1; %0.265 this is a corrective factor (obtained by try-and-error) that led to obtain a similar number of fragments wrt Somma
% phi=f_corrective*1/2*pi*sigma*v; % new (from Somma) the factor 1/2 is not correct
phi = f_corrective*pi*sigma*v; % new (from Somma)

%==========================================================================
 % Initial conditions
%==========================================================================

% Constant values
x0_S = ones(1,N_shell)*1e3 *1;
x0_D = ones(1,N_shell)*1e2 *1;
x0_N = ones(1,N_shell)*1e1 *1;
x0_lambda = ones(1,N_shell)*1e1 *0;
% x0_lambda = 10.^( linspace(3,1,length(R02(2:end))) ) *1;

% % Distributions
% mu_S=500;sig_S=300^2;S0=1000;
% mu_D=450;sig_D=300^2;D0=500;
% mu_N=300;sig_N=150^2;N0=400;
% mu_lam_S=350;sig_lam_S=300^2;lam0_S=3000;
% x0_lambda = lam0_S*exp(-(R02(2:end)-mu_lam_S).^2/sig_lam_S); 
% x0_S = S0*exp(-(R02(2:end)-mu_S).^2/sig_S);
% x0_D = D0*exp(-(R02(2:end)-mu_D).^2/sig_D);
% x0_N = N0*exp(-(R02(2:end)-mu_N).^2/sig_N);

%==========================================================================
% Propagation settings
%==========================================================================

tf_ss = 100; % total time
n_steps = 1000; % time steps
options_ode = odeset('reltol', 1e-2,'abstol', 1e-2);

%==========================================================================
% Model
%==========================================================================

S = sym('S',[N_shell,1]);
D = sym('D',[N_shell,1]);
N = sym('N',[N_shell,1]);
lam = sym('lam',[N_shell,1]);

F107_min = 70;
Ap_min = 2;
F107_max = 200;
Ap_max = 8;
F107_middle = 150;
Ap_middle = 5;

for k=1:N_shell 

    T = 900 + 2.5 *( F107_middle - 70 ) + 1.5* Ap_middle;
    m = 27 - 0.012* ( R02(k:k+1) - 200 ); %180 < h(km) < 1000 (is the model valid above 1000 km? if yes, use this, otherwise uncomment the lines below)
    H = T ./ m;
    
    if k<N_shell
        n0=(N(k+1));
        n_upper=D(k+1);
%         rho_k1= 6e-10* exp ( - ( R02(k+1) - 175 ) / H(2) ) ; % to include solar flux
        rho_k1 = densityexp(R02(k+1)); % no solar flux
        rvel_upper_D=-rho_k1*beta(2)*sqrt(mu*(R0(k+1)))*(24*3600*365.25);% drag flux
        rvel_upper_N=-rho_k1*beta(3)*sqrt(mu*(R0(k+1)))*(24*3600*365.25);%drag flux
        
    else
        % ASK MILES HERE!
        % n_upper = D(k);
        % n0 = N(k);
        n_upper = 0;
        n0 = 0;

        rho_k1 = densityexp((R02(k+1)+(R02(k+1)-R02(k)))); % no solar flux
        m = 27 - 0.012* ( (R02(k+1)+(R02(k+1)-R02(k))) - 200 ); %180 < h(km) < 1000 (is the model valid above 1000 km? if yes, use this, otherwise uncomment the lines below)
        H_top = T ./ m;

        rvel_upper_D=-rho_k1*beta(2)*sqrt(mu*(R0(k+1)+R0(k+1)-R0(k)))*(24*3600*365.25);% drag flux
        rvel_upper_N=-rho_k1*beta(3)*sqrt(mu*(R0(k+1)+R0(k+1)-R0(k)))*(24*3600*365.25);%drag flux
    end
%     rhok= (6e-10)* exp ( - ( R02(k) - 175 ) / H(1) ) ; % to include solar flux
    rhok = densityexp(R02(k)); % no solar flux
    rvel_current_D=-rhok*beta(2)*sqrt(mu*(R0(k)))*(24*3600*365.25);
    rvel_current_N=-rhok*beta(3)*sqrt(mu*(R0(k)))*(24*3600*365.25);

    phi_SS = phi(1,1)/V(k);
    phi_SD = phi(1,2)/V(k);
    phi_SN = phi(1,3)/V(k);
    phi_DD = phi(2,2)/V(k);
    phi_DN = phi(2,3)/V(k);
    phi_NN = phi(3,3)/V(k); 
    K0_SS = K0(1,1);
    K0_SD = K0(1,2);
    K0_SN = K0(1,3);
    K0_DD = K0(2,2);
    K0_DN = K0(2,3);
    K0_NN = K0(3,3);

    Eq1(k)=lam(k)-S(k)/Dt-(delta+alpha)*(phi_SN*N(k)*S(k)+phi_SD*D(k)*S(k))...
        -alpha_active*phi_SS*S(k)^2;

    Eq2(k)=(1-P)*S(k)/Dt+...
        delta*phi_SD*D(k)*S(k)+delta*phi_SN*N(k)*S(k)-(phi_DN*N(k)*D(k)+phi_DD*D(k)*D(k))...%creation
        +n_upper*rvel_upper_D/Dhu+D(k)*rvel_current_D/Dhl;%flux

    Eq3(k)=(K0_SN*phi_SN*N(k)*(alpha*S(k))+K0_SD*phi_SD*D(k)*(alpha*S(k)))+(K0_DN*phi_DN*N(k)*D(k)+...
        K0_DD*phi_DD*D(k)*D(k))+alpha_active*K0_SS*phi_SS*S(k)^2+phi_NN*K0_NN*N(k)^2 ...%creation 
            +n0*rvel_upper_N/Dhu+N(k)*rvel_current_N/Dhl; %flux [+n0/tau_res_disc(k+1,3)-N(k)/tau_res_disc(k,3)]
end

var=[S;D;N];
f3=[Eq1.';Eq2.';Eq3.'];
fun3= matlabFunction(f3,'Vars',{var,lam});
func_temp = @(x) fun3(x,x0_lambda');
func=@(t,x) func_temp(x);

%==========================================================================
% Propagation
%==========================================================================

x0 = [x0_S, x0_D, x0_N x0_lambda];

x0 = [    19
         204
         295
         382
           0
          19
         204
         295
         382
           0
          62
         443
         699
        1758
           0
           x0_lambda'];
tspan1 = linspace(0,tf_ss,n_steps);
[t_prop,x_prop] = ode15s(func, tspan1, x0(1:3*N_shell), options_ode);

S_prop = x_prop(:,1:N_shell);
D_prop = x_prop(:,N_shell+1:2*N_shell);
N_prop = x_prop(:,2*N_shell+1:3*N_shell);
N_tot = S_prop + D_prop + N_prop;
N_tot_sum = sum(N_tot,2);
S_sum = sum(S_prop,2);
D_sum = sum(D_prop,2);
N_sum = sum(N_prop,2);

vpa(f3,4)

%==========================================================================
% Plots
%==========================================================================

figure;
hold on; grid on;
plot(R02(2:end),x0(1:N_shell),'LineWidth',2);
plot(R02(2:end),x0(N_shell+1:2*N_shell),'LineWidth',2);
plot(R02(2:end),x0(2*N_shell+1:3*N_shell),'LineWidth',2);
plot(R02(2:end),x0(3*N_shell+1:4*N_shell),'LineWidth',2);
xlabel("Altitude [km]");ylabel('Number')
legend("S_{0}","D_{0}","N_{0}","\lambda_{0}")
title("I.C.")
set(gca, 'FontSize', 16)

figure;
hold on; grid on;
plot(t_prop,N_tot_sum,'LineWidth',2)
plot(t_prop,S_sum,'LineWidth',2);
plot(t_prop,D_sum,'LineWidth',2);
plot(t_prop,N_sum,'LineWidth',2);
xlabel("Time [years]");ylabel('Number of objects')
legend("total","S","D","N")
title("Population VS Time")
set(gca, 'FontSize', 16)

figure;
subplot(3,1,1)
surf(t_prop', R02(2:end), transpose(S_prop(:,:,:)),'edgecolor','none'); % Note the transpose of the data matrix
c = colorbar;
c.Label.String = 'Number of Objects';
clim([0 max(S_prop(:))])
zlabel('Number of Objects')
xlabel('Time (years)')
ylabel('Altitude (km)')
xlim([0 max(t_prop)])
ylim([min(R02(2:end)) max(R02(2:end))])
grid on
view (0, 90)
title("S")
set(gca, 'FontSize', 16)

% figure;
subplot(3,1,2)
surf(t_prop', R02(2:end), transpose(D_prop),'edgecolor','none'); % Note the transpose of the data matrix
c = colorbar;
c.Label.String = 'Number of Objects';
clim([0 max(D_prop(:))])
zlabel('Number of Objects')
xlabel('Time (years)')
ylabel('Altitude (km)')
xlim([0 max(t_prop)])
ylim([min(R02(2:end)) max(R02(2:end))])
grid on
view (0, 90)
title("D")
set(gca, 'FontSize', 16)

% figure;
subplot(3,1,3)
surf(t_prop', R02(2:end), transpose(N_prop),'edgecolor','none'); % Note the transpose of the data matrix
c = colorbar;
c.Label.String = 'Number of Objects';
clim([0 max(N_prop(:))])
zlabel('Number of Objects')
xlabel('Time (years)')
ylabel('Altitude (km)')
xlim([0 max(t_prop)])
ylim([min(R02(2:end)) max(R02(2:end))])
grid on
view (0, 90)
title("N")
set(gca, 'FontSize', 16)
