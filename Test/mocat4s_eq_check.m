function [f4, OUT] = mocat4s_eq_check()

% close all
% clear 
% clc

addpath(genpath(pwd)) % add subfolders to path
% rmpath(genpath(pwd)) % remove subfolders to path

Si = [17 373 1232 233 1];
Sui = [10 54 2439 27 1];
Di = [19 204 295 382 1];
Ni = [62 443 699 1758 1];

VAR = MOCAT4S_VAR_Cons();

% Integrate over 200 years with constant yearly launch rate
tf = 100; %Years

x0 = [Si';Di';Ni';Sui'];

lam(:,1) = .05 * Si' * 0;
lam(:,2) = .05 * Sui' * 0;

VAR.N_step = 1000;
tspan = linspace(0,tf,VAR.N_step);
tic
OUT = MOCAT4S(tspan, x0, lam, VAR);
toc

X = [OUT.S'; OUT.D'; OUT.N'; OUT.Su']';
tf = tspan';
LAM1 = ones(size(tf,1),VAR.N_shell).*lam(:,1)';
LAM2 = ones(size(tf,1),VAR.N_shell).*lam(:,2)';

% %% Integrate over multiple time periods with feedback controller yearly determined launch rate
% VAR.N_step = 2;
% X = x0';
% LAM1 = lam(:,1)' * 0;
% LAM2 = lam(:,2)' * 0;
% tf = [0:5:tf]; %Years
% for t = 1:length(tf)-1
%     disp(tf(t))
%     tspan = linspace(tf(t),tf(t+1),VAR.N_step);
%     OUT = MOCAT4S(tspan, x0, lam, VAR);
%     % Launch rates for NEXT time step based on states in previous time step
%     lam(:,1) = 0.05 * OUT.S(end,:)';
%     lam(:,2) = 0.05 * OUT.Su(end,:)';
%     % lam(:,1) = .05 * Si'; % same as line 43-44
%     % lam(:,2) = .05 * Sui';
%     x0 = [OUT.S(end,:)'; OUT.D(end,:)'; OUT.N(end,:)'; OUT.Su(end,:)'];
%     X = [X;x0'];
%     LAM1 = [LAM1;lam(:,1)'];
%     LAM2 = [LAM2;lam(:,2)'];
% end

%% Graph outputs of iterative model run
S = X(:,1:VAR.N_shell);
D = X(:,VAR.N_shell+1:2*VAR.N_shell);
N = X(:,2*VAR.N_shell+1:3*VAR.N_shell);
Su = X(:,3*VAR.N_shell+1:4*VAR.N_shell);
N_tot = S+D+N+Su; % total population at each time instant for each shell
% 
% figure
% subplot(1,3,1)
% hold on
% grid on
% set(gca,'FontSize',14)
% plot(VAR.R02(2:end),S(1,:),'k*-','LineWidth',2)
% plot(VAR.R02(2:end),S(end,:),'*-','LineWidth',2)
% xlabel('Altitude [km]')
% ylabel('Slotted satellites')
% legend('S_0','S_f')
% 
% subplot(1,3,2)
% hold on
% grid on
% set(gca,'FontSize',14)
% plot(VAR.R02(2:end),D(1,:),'k*-','LineWidth',2)
% plot(VAR.R02(2:end),D(end,:),'*-','LineWidth',2)
% xlabel('Altitude [km]')
% ylabel('Derelicts')
% legend('D_0','D_f')
% 
% subplot(1,3,3)
% hold on
% grid on
% set(gca,'FontSize',14)
% plot(VAR.R02(2:end),N(1,:),'k*-','LineWidth',2)
% plot(VAR.R02(2:end),N(end,:),'*-','LineWidth',2)
% xlabel('Altitude [km]')
% ylabel('Debris')
% legend('N_0','N_f')
% 
% figure
% semilogy(VAR.R02(2:end),S(end,:),'LineWidth',2)
% hold on
% semilogy(VAR.R02(2:end),Su(end,:),'LineWidth',2)
% grid on
% set(gca,'FontSize',14)
% xlabel('Altitude [km]')
% ylabel('Number of objects')
% legend('S','Su') %
% title('Final time')
% 
% figure
% hold on
% grid on
% set(gca,'FontSize',14)
% xlabel('Altitude [km]')
% ylabel('Launch Rate S (Number of objects/year)')
% for t = 1:length(tf)-1
%     plot(VAR.R02(2:end), LAM1(t,:))
% end
% 
% figure
% hold on
% grid on
% set(gca,'FontSize',14)
% xlabel('Altitude [km]')
% ylabel('Launch Rate Su (Number of objects/year)')
% for t = 1:length(tf)-1
%     plot(VAR.R02(2:end), LAM2(t,:))
% end
% 
% figure 
% hold on
% grid on
% set(gca,'FontSize',14)
% plot(tf,sum(N_tot,2),'-','LineWidth',2)
% plot(tf,sum(S,2),'-','LineWidth',2)
% plot(tf,sum(D,2),'-','LineWidth',2)
% plot(tf,sum(N,2),'-','LineWidth',2)
% plot(tf,sum(Su,2),'-','LineWidth',2)
% xlabel('Time [years]')
% ylabel('Number of objects')
% legend('total','S','D','N','Su')
% title('Population vs time')
% 
% 
% x = tf';
% y = VAR.R02(2:end);
% [Xi,Yi] = meshgrid(x,y);
% 
% figure
% subplot(4,1,1)
% hold on
% % grid on
% set(gca,'FontSize',14)
% % surf(Xi,Yi,S','edgecolor','none')
% surf(tf', VAR.R02(2:end), transpose(S),'edgecolor','none'); 
% c = colorbar;
% c.Label.String = 'S';
% xlabel('Time [years]')
% ylabel('Altitude [km]')
% ylim([VAR.R02(2) VAR.R02(end)])
% 
% subplot(4,1,2)
% hold on
% % grid on
% set(gca,'FontSize',14)
% surf(Xi,Yi,Su','edgecolor','none')
% c = colorbar;
% c.Label.String = 'S_u';
% xlabel('Time [years]')
% ylabel('Altitude [km]')
% ylim([VAR.R02(2) VAR.R02(end)])
% 
% subplot(4,1,3)
% hold on
% % grid on
% set(gca,'FontSize',14)
% surf(Xi,Yi,D','edgecolor','none')
% c = colorbar;
% c.Label.String = 'D';
% xlabel('Time [years]')
% ylabel('Altitude [km]')
% ylim([VAR.R02(2) VAR.R02(end)])
% 
% subplot(4,1,4)
% hold on
% % grid on
% set(gca,'FontSize',14)
% surf(Xi,Yi,N','edgecolor','none')
% c = colorbar;
% c.Label.String = 'N';
% xlabel('Time [years]')
% ylabel('Altitude [km]')
% ylim([VAR.R02(2) VAR.R02(end)])


%% FUNCTIONS
function [VAR] = MOCAT4S_VAR_Cons()

% Constants
rad = pi/180;
mu_km=398600.4415;%km^3/s^2
mu=3.986004418e14;%meters^3/s^2
re = 6378.1366; % [km]
years_s = 365*24*3600;

% Parameters needed for all MOCAT-SSEMs
VAR.N_shell = 5;
VAR.h_max = 900;
VAR.h_min = 200;
VAR.N_step = 1000; % how many steps to calculate between t0 and tf
N_shell = VAR.N_shell;
h_max = VAR.h_max;
h_min = VAR.h_min;
R0=linspace(h_min,h_max,N_shell+1);
R02=linspace(h_min,h_max,N_shell+1);
VAR.deltaH=R02(2)-R02(1); % thickness of the shell [km]

% Parameters needed for MOCAT-4S
VAR.alpha = 2e-3; % fraction of collisions that satellites fail to avoid.
VAR.alpha_active= 1e-3; % fraction of collisions occuring between active satellites 
VAR.Dt = 5; % operation lifetime: after those years, active satellites become derelict [years]
VAR.slot_par = 1; % slotting effectiveness parameter in [0,1] where 1 = no cols. between slotted satellites.
VAR.delta = 10; %10
VAR.P = 0.9;
VAR.Cd = 2.2; % Coefficent of drag for drag purposes
VAR.mass      = [223 223 0.640 223];      % median mass values from MASTER [kg] SOMMA
VAR.A         = [1.741 1.741 0.020 1.741]; % median area values from MASTER [m^2] SOMMA
VAR.diameter  = [1.490 1.490 0.180 1.490]; % median diameter values from MASTER [m] SOMMA
VAR.v_imp = 10; % impact velocity [km/s]
VAR.LC = 0.1; % minimum size of fragments [m]
VAR.deg_intrinsic = 45; 

VAR.sep_dist_method = "distance";% either 'distance' or 'angle'. Sets minimum separation distance for intrinsic capacity constraint.
VAR.sep_angle = 0.1; % minimum angular separation distance [deg]. Only used if sep_dist_method VAR.sep_dist_method = 'angle'
VAR.sep_dist = 5; % minimum separation distance [km]. Only used if sep_dist_method VAR.sep_dist_method = 'distance'

% VAR.input_pop = 'TLE';  % 'TLE' or 'distribution'.  TLE pulls binned present day objects. 'distribution' uses an arbitary hardcoded population distribution.
% Generate TLE binned appropriately using the getTLEBins python method
% TODO: 
% savedir = py.tletobins.getTLEBins(VAR.deltaH, VAR.h_min, VAR.h_max, savedir = None, graph = False)
% VAR.filename_N = fullfile(savedir, 'Counts_DEBRIS_bins_' +num2str(VAR.deltaH) + '.csv';
% VAR.filename_S = fullfile(savedir, 'Counts_PAYLOADslot_bins_' +num2str(VAR.deltaH) + '.csv'
% VAR.filename_Su = fullfile(savedir, 'Counts_PAYLOADunslot_bins_' +num2str(VAR.deltaH) + '.csv'
% VAR.filename_D = fullfile(savedir, 'Counts_DERELICT_bins_' +num2str(VAR.deltaH) + '.csv'

% Change this.
% VAR.filename_N = '..\MOCAT-4S\x0_TLE\Counts_DEBRIS_bins_35.csv';
% VAR.filename_S = '..\MOCAT-4S\x0_TLE\Counts_PAYLOADslot_bins_35.csv';
% VAR.filename_Su = '..\MOCAT-4S\x0_TLE\Counts_PAYLOADunslot_bins_35.csv';
% VAR.filename_D = '..\MOCAT-4S\x0_TLE\Counts_DERELICT_bins_35.csv';
% 
% VAR.filename_N = '../MOCAT-4S/x0_TLE/Counts_DEBRIS_bins_35.csv';
% VAR.filename_S = '../MOCAT-4S/x0_TLE/Counts_PAYLOADslot_bins_35.csv';
% VAR.filename_Su = '../MOCAT-4S/x0_TLE/Counts_PAYLOADunslot_bins_35.csv';
% VAR.filename_D = '../MOCAT-4S/x0_TLE/Counts_DERELICT_bins_35.csv';


R1=R0;
R0=(re+R0)*1000;
% V=4*pi*R0.^2*deltaH*1000;
VAR.V=4/3*pi*(diff(R0.^3)); % volume of the shells [m^3]
VAR.v=VAR.v_imp*1000*(24*3600*365.25);% impact velocity [m/year]
VAR.Dhl=VAR.deltaH*1000;
VAR.Dhu=-VAR.deltaH*1000;
Cd=VAR.Cd;

% physical characteristics
n_obj = length(VAR.mass);
radii = VAR.diameter/2;
VAR.area_mass = VAR.A./VAR.mass; % area/mass ratio  [# shells, # object types]
beta = VAR.Cd*VAR.area_mass; % ballistic coefficient [# shells, # object types]

% NASA standard breakup models implementation for collisions
v_imp = VAR.v_imp; % impact velocity [km/s]
LC = VAR.LC; % minimum size of fragments [m]
n_f_catastrophic = @(M1,M2) 0.1*LC^(-1.71)*(M1+M2)^(0.75); % number of fragments generated during a catastrophic collision (NASA standard break-up model). M is the sum of the mass of the objects colliding in kg
n_f_damaging = @(M1,M2) 0.1*LC^(-1.71)*(min(M1,M2)*v_imp^2)^(0.75); % number of fragments generated during a non-catastrophic collision (improved NASA standard break-up model: takes into account the kinetic energy). M is the mass of the less massive object colliding in kg
K0 = zeros(n_obj,n_obj); % number of fragments generated in each collisions among the species for each altitude bin
ind_obj_catastrophic = 3; % hypothesis of catastrophic collisions only between satellites and derelict satellites
sigma = zeros(n_obj,n_obj);
for q = 1:length(radii)
    mass_q = VAR.mass(q);
    for qq = 1:length(radii)
        sigma(q,qq) = (radii(q)+radii(qq))^2;
        mass_qq = VAR.mass(qq);
        if q~=ind_obj_catastrophic && qq~=ind_obj_catastrophic
            K0(q,qq) = n_f_catastrophic(mass_q,mass_qq); % catastrophic collision
        else
            K0(q,qq) = n_f_damaging(mass_q,mass_qq); % (non-catastrophic) damaging collision
        end
    end
end

% phi=sigma*v; % adopted before 
f_corrective = 1; %0.265 this is a corrective factor (obtained by try-and-error) that led to obtain a similar number of fragments wrt Somma
% phi=f_corrective*1/2*pi*sigma*v; % new (from Somma) the factor 1/2 is not correct
phi=f_corrective*pi*sigma*VAR.v; % new (from Somma)

%% initial conditions

%input_pop = 'TLE';
% input_pop = 'distribution';

% if strcmp(VAR.input_pop,'TLE')
%     %%% MODIFY FOR MORE SHELLS
%     % load initial conditions from TLEs
% 
%     filename_N = VAR.filename_N;
%     filename_S = VAR.filename_S;
%     filename_Su = VAR.filename_Su;
%     filename_D = VAR.filename_D;
% 
%     Ni = readmatrix(filename_N);
%     Ni = Ni(:,2)';
% 
%     Si = readmatrix(filename_S);
% %   Si = Si(:,3)';
%     Si = Si(:,2)';
% 
%     Sui = readmatrix(filename_Su);
%     Sui = Sui(:,2)';
% 
%     Di = readmatrix(filename_D);
%     Di = Di(:,2)';
% %   Di = floor(0.5*(Si+Sui));
% 
% elseif strcmp(VAR.input_pop,'distribution')
% 
%     % previous initial conditions
% 
%     mu_S=500;sig_S=300^2;S0=1000;
%     mu_Su=650;sig_Su=300^2;Su0=850;
%     mu_D=450;sig_D=300^2;D0=500;
%     mu_N=300;sig_N=150^2;N0=400;
%     mu_lam_S=350;sig_lam_S=300^2;lam0_S=1000; %lam0=3000
%     mu_lam_Su=350;sig_lam_Su=300^2;lam0_Su=1000; %lam0=3000
%     lam_model_S=lam0_S*exp(-(R02(2:end)-mu_lam_S).^2/sig_lam_S); 
%     lam_model_Su=lam0_Su*exp(-(R02(2:end)-mu_lam_Su).^2/sig_lam_Su); 
% 
%     Si=S0*exp(-(R02(2:end)-mu_S).^2/sig_S);
%     Sui=Su0*exp(-(R02(2:end)-mu_Su).^2/sig_Su);
%     Di=D0*exp(-(R02(2:end)-mu_D).^2/sig_D);
%     Ni=N0*exp(-(R02(2:end)-mu_N).^2/sig_N);
% 
% end

% start computation

options = odeset('reltol', 1.e-4,'abstol', 1.e-4);

VAR.x0 = [Si';Di';Ni';Sui'];
VAR.options = options;
%VAR.h = h;
VAR.beta = beta;
VAR.mu = mu;
VAR.K0 = K0;
VAR.phi = phi;
VAR.R0 = R0;
VAR.R02 = R02;

launch_rate = 'generic';
VAR.S_Su = 0.1; %ratio of unslotted to slotted spacecraft with generic or distribution model
MOCAT = '4S'; 
VAR.launch_rate = launch_rate;

end


function [OUT] = MOCAT4S(tspan, x0, lam, VAR)
% tspan = times for integration (in years)
% x0 = current population VAR.N_shell x 4 in order S, D, N,  Su
% lam matrix S, Su per bin so VAR.N_shell x 2 in order S Su
%
warning('off')

%% propagation

S = sym('S_',[VAR.N_shell,1]);
D = sym('D_',[VAR.N_shell,1]);
N = sym('N_',[VAR.N_shell,1]);
Su = sym('Su_',[VAR.N_shell,1]);
lam_S = sym('lam_S',[VAR.N_shell,1]);
lam_Su = sym('lam_Su',[VAR.N_shell,1]);

lam1_S = lam(:,1);
lam1_Su = lam(:,2);

% compute derivatives

for k=1:VAR.N_shell

    if k<VAR.N_shell

        n0=(N(k+1));
        n_upper=D(k+1);

        rho_k1 = densityexp(VAR.R02(k+1)); % no solar flux
        rvel_upper_D=-rho_k1*VAR.beta(2)*sqrt(VAR.mu*(VAR.R0(k+1)))*(24*3600*365.25);% drag flux
        rvel_upper_N=-rho_k1*VAR.beta(3)*sqrt(VAR.mu*(VAR.R0(k+1)))*(24*3600*365.25);%drag flux

    else
        % ASSUMPTION: No flux coming down from highest shell.
        n_upper = 0; 
        n0 = 0;
        rho_k1 = densityexp((VAR.R02(k+1)+(VAR.R02(k+1)-VAR.R02(k)))); % no solar flux
        rvel_upper_D=-rho_k1*VAR.beta(2)*sqrt(VAR.mu*(VAR.R0(k+1)+VAR.R0(k+1)-VAR.R0(k)))*(24*3600*365.25);% drag flux
        rvel_upper_N=-rho_k1*VAR.beta(3)*sqrt(VAR.mu*(VAR.R0(k+1)+VAR.R0(k+1)-VAR.R0(k)))*(24*3600*365.25);%drag flux

    end

    rhok = densityexp(VAR.R02(k)); % no solar flux
    rvel_current_D=-rhok*VAR.beta(2)*sqrt(VAR.mu*(VAR.R0(k)))*(24*3600*365.25);
    rvel_current_N=-rhok*VAR.beta(3)*sqrt(VAR.mu*(VAR.R0(k)))*(24*3600*365.25);

    % use different probability of collision and generated number of
    % fragments according to the species colliding

    phi_SS = VAR.phi(1,1)/VAR.V(k);
    phi_SD = VAR.phi(1,2)/VAR.V(k);
    phi_SN = VAR.phi(1,3)/VAR.V(k);
    phi_SSu = VAR.phi(1,4)/VAR.V(k);
    phi_DD = VAR.phi(2,2)/VAR.V(k);
    phi_DN = VAR.phi(2,3)/VAR.V(k);
    phi_DSu = VAR.phi(2,4)/VAR.V(k);
    phi_NN = VAR.phi(3,3)/VAR.V(k);
    phi_NSu = VAR.phi(3,4)/VAR.V(k);
    phi_SuSu = VAR.phi(4,4)/VAR.V(k);
    
    K0_SS = VAR.K0(1,1);
    K0_SD = VAR.K0(1,2);
    K0_SN = VAR.K0(1,3);
    K0_SSu = VAR.K0(1,4);
    K0_DD = VAR.K0(2,2);
    K0_DN = VAR.K0(2,3);
    K0_DSu = VAR.K0(2,4);
    K0_NN = VAR.K0(3,3);
    K0_NSu = VAR.K0(3,4);
    K0_SuSu = VAR.K0(4,4);

    A0 = lam_S(k);
    A1 = -1/VAR.Dt;
    A2 = -(VAR.delta+VAR.alpha)*phi_SN;
    A3 = -(VAR.delta+VAR.alpha)*phi_SD;
    A4 = -(1-VAR.slot_par)*VAR.alpha_active*phi_SS;
    A5 = -VAR.alpha_active*phi_SSu*S(k)/(Su(k)+S(k)); 

    Eq1(k) = A0+A1*S(k)+A2*S(k)*N(k)+A3*S(k)*D(k)+A4*S(k)^2 +A5*S(k)*Su(k);

    B0 = +n_upper*rvel_upper_D/VAR.Dhu;
    B1 = (1-VAR.P)/VAR.Dt;%S(k)
    B2 = rvel_current_D/VAR.Dhl;%D(k)
    B3 = +VAR.delta*phi_SD;%D(k)*S(k)
    B4 = -phi_DN;%(N(k))*D(k)
    B5 = -phi_DD;%(D(k))*D(k)
    B6 = +VAR.delta*phi_SN;%N(k)*S(k)
    B7 = (1-VAR.P)/VAR.Dt;%Su(k)
    B8 = +VAR.delta*phi_DSu;%D(k)*Su(k)
    B9 = +VAR.delta*phi_NSu;%N(k)*Su(k)

    Eq2(k) = B1*S(k)+B3*D(k)*S(k)+B6*S(k)*N(k)+B4*D(k)*N(k)+B5*D(k)^2+B0+B2*D(k) +B7*Su(k)+B8*D(k)*Su(k)+B9*N(k)*Su(k);

    C0 = +n0*rvel_upper_N/VAR.Dhu;
    C1 = rvel_current_N/VAR.Dhl;%N(k)
    C2 = +phi_SD*VAR.alpha*K0_SD;%D(k)*S(k) %new
    C3 = K0_DN*phi_DN;%(N(k))*D(k)
    C4 = K0_DD*phi_DD;%(D(k))*D(k)
    C5 = phi_SN*VAR.alpha*K0_SN;%(N(k))*S(k)
    C6 = (1-VAR.slot_par)*VAR.alpha_active*phi_SS*K0_SS;%S(k)^2
    C7 = phi_NN*K0_NN;%N(k)^2 %new
    C8 = +phi_DSu*VAR.alpha*K0_DSu;%D(k)*Su(k) %new
    C9 = phi_NSu*VAR.alpha*K0_NSu;%(N(k))*S(k)
    C10 = VAR.alpha_active*phi_SuSu*K0_SuSu;
    C11 = +VAR.alpha_active*K0_SSu*phi_SSu; % S(k)*Su(k) 

    Eq3(k) = C5*N(k)*S(k)+C2*D(k)*S(k)+C3*N(k)*D(k)+C4*D(k)^2+C6*S(k)^2+C7*N(k)^2+C0+C1*N(k) +C8*D(k)*Su(k)+C9*N(k)*Su(k)+C10*Su(k)^2+C11*S(k)*Su(k);

    A0_Su = lam_Su(k);
    A1_Su = -1/VAR.Dt;
    A2_Su = -(VAR.delta+VAR.alpha)*phi_NSu;
    A3_Su = -(VAR.delta+VAR.alpha)*phi_DSu;
    A4_Su = -VAR.alpha_active*phi_SuSu;
    A5_Su = -VAR.alpha_active*phi_SSu*Su(k)/(Su(k)+S(k)); 

    Eq4(k) = A0_Su+A1_Su*Su(k)+A2_Su*Su(k)*N(k)+A3_Su*Su(k)*D(k)+A4_Su*Su(k)^2 +A5_Su*S(k)*Su(k);

end

var=[S;D;N;Su];
f4=[Eq1.';Eq2.';Eq3.';Eq4.'];
f4 = subs(f4, [lam_S, lam_Su], [lam1_S,lam1_Su]);
fun4= matlabFunction(f4,'Vars',{var});
func_temp = @(x) (fun4(x));
func=@(t,x) func_temp(x);
options_ode = odeset('reltol',1.e-4,'abstol', 1.e-4);
[T,X] = ode15s(func, tspan, x0, options_ode);

disp("Equations-4S")
vpa(f4,4)

% tspan
% T
% X

%%%
S = X(:,1:VAR.N_shell);
D = X(:,VAR.N_shell+1:2*VAR.N_shell);
N = X(:,2*VAR.N_shell+1:3*VAR.N_shell);
Su = X(:,3*VAR.N_shell+1:4*VAR.N_shell);

% N_tot = S+D+N+Su; % total population at each time instant for each shell
% dens_shell = N_tot./(VAR.V*1e-9); % total spatial density over time along all the shells [obj/km^3]
% dens_tot_in = N_tot(1,:)./(VAR.V*1e-9); % Initial total spatial density along all the shells [obj/km^3]
% dens_tot_end = N_tot(end,:)./(VAR.V*1e-9); % Final total spatial density along all the shells [obj/km^3]
% N_tot_real_S = sum(S(end,:)); % total number of slotted satellites 
% N_tot_real_Su = sum(Su(end,:)); % total number of unslotted satellites 
% N_tot_real = N_tot_real_S+N_tot_real_Su; % total number of satellites 
% N_diff_debris = (sum(N(end,:))-sum(N(1,:)))/sum(N(1,:)); % difference in the number of debris
% dens_diff_debris = (N(end,:)-N(1,:))./(VAR.V*1e-9); % difference in the number of debris [N/km^3]
% N_der = diff(N(end-1:end,:))/(T(end)-T(end-1));

OUT = struct;
OUT.S = S;
OUT.D = D;
OUT.N = N;
OUT.Su = Su;
OUT.T = T;

end
end