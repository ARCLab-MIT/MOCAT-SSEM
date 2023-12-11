clearvars; close all; clc

%% Plot settings
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',22);
set(0,'defaultLineLineWidth',1.8);
set(0,'defaultLineMarkerSize',7);
rgb_c = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],...
         [0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840]};


%% Read MASTER tables
% [alt [km], diametere, collision fragments, explosion fragments, LMRO, total]
T_LMRO = table2cell(readtable('IPOP_200km_2000km.xlsx','Sheet','LMRO','Range','A2:F36180'));    
T_N = table2cell(readtable('IPOP_200km_2000km.xlsx','Sheet','N','Range','A2:E181'));  
T = {T_LMRO, T_N};


% LMRO:
% diameter in basic setting     0.1 - 100   [m]
% diameter in 2D spectrum       0.1 - 100   classes -200 (log flagged)
% altitude in 2D spectrum       200 - 2000  classes +10
% 
% COLL+EXPL:
% diameter in basic setting     0.1 - 1     [m]
% diameter in 2D spectrum       0.1 - 1     classes -200 (log flagged)
% altitude in 2D spectrum       200 - 2000  classes +10


%% Create matrices

% Altitudes
h_min = 200;        % [km]
h_max = 2000;       % [km]
N_shells_master = 181;     % for 10km bin size
R_master = linspace(h_min,h_max,N_shells_master);
alt_meas = 205:10:1995;
re = 6378.1;
R0 = (re+R_master)*1000;
V = 4/3*pi*(diff(R0.^3));   % volume of the shells [m^3]
Vkm = V*1e-9;               % volume of the shells [km^3]

% Number of diameters from MASTER discretisation
N_diam = 200;

% Vector of diameter
d_vect = transpose([T_LMRO{1:N_diam,2}]);

% Create matrices
pop_lmro = zeros(N_shells_master-1,N_diam+1);   % add one column for the altitude
pop_n = zeros(N_shells_master-1,2);      % add one column for the altitude
k = 1;

for i = 1:N_shells_master-1
    for j = 2:N_diam+1
        pop_lmro(i,j) = [T_LMRO{k,end}];
%         pop_n(i,j) = [T_N{k,end}];
        k = k+1;
    end
    k = k+1;
    pop_n(i,2) = [T_N{i,end}];
end

% From spatial density to number of objects
pop_lmro = pop_lmro .* Vkm';
pop_n = pop_n .* Vkm';

% Add the altitude column corresponding to MASTER discretisation
pop_lmro(:,1) = alt_meas;
pop_n(:,1) = alt_meas;


%% PLOT
figure();
% subplot(1,2,1)
[X_surf,Y_surf] = meshgrid(d_vect,alt_meas);
h = surf(X_surf,Y_surf,pop_lmro(:,2:end));
xline(2,'r',{'$d_{\rm{cutoff}}$'},'interprete','latex','LineWidth',2.5,...
    'labelorientation','horizontal','fontsize',16)
c = colorbar;
c.Title.Interpreter = 'latex';
c.Title.String = '[\# obj]';
c.Label.FontSize = 12;
title('LMRO')
set(h,'LineStyle','none')
ylabel('Altitude [km]'); xlabel('Diameter')
view(0,90)
ylim([200 2000])
set(gca, 'Xscale','log')

% subplot(1,2,2)
% h = surf(X_surf,Y_surf,pop_n(:,2:end));
% c = colorbar;
% c.Title.Interpreter = 'latex';
% c.Title.String = '[\# obj]';
% c.Label.FontSize = 12;
% title('Debris')
% set(h,'LineStyle','none')
% view(0,90)
% ylabel('Altitude [km]'); xlabel('Diameter')
% ylim([200 2000])
% % xlim([2 20]); axis([0 1e-11])                %%%%%%%% there are some fragments with RB dimensions
% set(gca, 'Xscale','log')

%% Construct the initial population

% Cutoff diameter
d_cutoff = 2;       % [m]

% Parameters
N_species = 4;
ind_S = 2;  ind_D = 3;  ind_N = 4;  ind_B = 5;
bin_pop = zeros(N_shells_master-1, N_species+1);
bin_pop(:,1) = transpose(205:10:1995);

% Subdivide the species
for i = 1:N_shells_master-1
    for j = 1:N_diam
        if d_vect(j) <= d_cutoff
            bin_pop(i,ind_S) = bin_pop(i,ind_S) + pop_lmro(i,j+1);
%             bin_pop(i,ind_N) = bin_pop(i,ind_N) + pop_n(i,j+1);                               % to consider rocket bodies also the big fragments
        else
            bin_pop(i,ind_B) = bin_pop(i,ind_B) + pop_lmro(i,j+1);
%             bin_pop(i,ind_S) = bin_pop(i,ind_S) + pop_lmro(i,j+1) + pop_n(i,j+1);             % to consider rocket bodies also the big fragments

        end
%         bin_pop(i,ind_N) = bin_pop(i,ind_N) + pop_n(i,j+1);                                     % to eliminate if we consider as B the big fragments
    end
end
bin_pop(:,ind_N) = pop_n(:,end);

% Distinguish between S and D
D_perc = 0.8;
bin_pop = round(bin_pop);
bin_pop(:,ind_D) = bin_pop(:,ind_S) * D_perc;
bin_pop(:,ind_S) = bin_pop(:,ind_S) - bin_pop(:,ind_D); 

%% PLOT, 10km alt bin
fig = figure;
subplot(2,2,1)
hold on; grid on; grid minor
plot(bin_pop(:,1),bin_pop(:,ind_S),'-k*')
ylabel('Active satellites'); xlabel('Altitude [km]'); xlim([h_min h_max])

subplot(2,2,2)
hold on; grid on; grid minor
plot(bin_pop(:,1),bin_pop(:,ind_D),'-k*')
ylabel('Derelict'); xlabel('Altitude [km]'); xlim([h_min h_max])

subplot(2,2,3)
hold on; grid on; grid minor
plot(bin_pop(:,1),bin_pop(:,ind_N),'-k*')
ylabel('Debris'); xlabel('Altitude [km]'); xlim([h_min h_max])

subplot(2,2,4)
hold on; grid on; grid minor
plot(bin_pop(:,1),bin_pop(:,ind_B),'-k*')
ylabel('Rocket bodies'); xlabel('Altitude [km]'); xlim([h_min h_max])

han=axes(fig,'visible','off');
han.Title.Visible='on';
title(han,['Initial population, h$_{shell} = $', num2str(R_master(2)-R_master(1))])


%% Discretisation, 50 km bin size
N_shells = 37;      % for 50km bin size
% N_shells = 73;      % for 25km bin size
R0 = linspace(h_min, h_max, N_shells);
% alt = flip(alt_meas);
IPop = zeros(N_shells, N_species+1);
IPop(:,1) = R0;
for k = 2:N_species+1
    j = 1;
    for i = 1:size(bin_pop,1)
        if alt_meas(i) < R0(j)
            IPop(j,k) = IPop(j,k) + bin_pop(i,k);
        else
            j = j + 1;
        end
    end
end

% Initial population 50km bin altitude
x0 = IPop(:,2:end);
x0 = x0(:);


%% PLOT, 50km alt bin
fig = figure;
subplot(2,2,1)
hold on; grid on; grid minor
plot(R0,IPop(:,ind_S),'-k*')
ylabel('Active satellites'); xlabel('Altitude [km]'); xlim([h_min h_max])

subplot(2,2,2)
hold on; grid on; grid minor
plot(R0,IPop(:,ind_D),'-k*')
ylabel('Derelict'); xlabel('Altitude [km]'); xlim([h_min h_max])

subplot(2,2,3)
hold on; grid on; grid minor
plot(R0,IPop(:,ind_N),'-k*')
ylabel('Debris'); xlabel('Altitude [km]'); xlim([h_min h_max])

subplot(2,2,4)
hold on; grid on; grid minor
plot(R0,IPop(:,ind_B),'-k*')
ylabel('Rocket bodies'); xlabel('Altitude [km]'); xlim([h_min h_max])

han=axes(fig,'visible','off');
han.Title.Visible='on'; 
title(han,['Initial population, h$_{shell} = $', num2str(R0(2)-R0(1))])



%% TEST
clc
close all

h_max = 2000;               % [km]
h_bin = 50;                 % [km]
dB_cutoff = 2;              % [m]
dN_cutoff = [0.05 0.08];    % [m]
species = '5N';             % define the species to consider: 3,4B,5N,(5S)
[IPop, IPop_interp] = BuiltIPop(h_max,h_bin,dB_cutoff,dN_cutoff,species);

% Initial vector
x0 = IPop(:,3:end);
x0 = x0(:);


% PLOT
ind_S = 3; ind_D = 4; ind_N = 5; ind_B = 6; ind_NU = 7;
fig = figure;
tl = tiledlayout(2,2);
nexttile
hold on; grid on; grid minor
plot(IPop(:,2),IPop(:,ind_S),'-k*')
plot(IPop_interp(:,1),IPop_interp(:,ind_S),'-r',LineWidth=1.2)
% ylabel('Active satellites [obj]'); xlabel('Altitude [km]'); 
title('Active satellites')
xlim([IPop(1,2) IPop(end,2)])

nexttile
hold on; grid on; grid minor
plot(IPop(:,2),IPop(:,ind_D),'-k*')
plot(IPop_interp(:,1),IPop_interp(:,ind_D),'-r',LineWidth=1.2)
% ylabel('Derelict [obj]'); xlabel('Altitude [km]'); 
xlim([IPop(1,2) IPop(end,2)])
title('Derelict')

nexttile
hold on; grid on; grid minor
plot(IPop(:,2),IPop(:,ind_N),'-k*')
plot(IPop_interp(:,1),IPop_interp(:,ind_N),'-r',LineWidth=1.2)
% ylabel('Debris [obj]'); xlabel('Altitude [km]'); 
title('Debris')
xlim([IPop(1,2) IPop(end,2)])
if strcmp(species,'5N')
    plot(IPop(:,2),IPop(:,end),'-k*')
    plot(IPop_interp(:,2),IPop_interp(:,end),'-b',LineWidth=1.2)
    legend({'','Tracked','',['Untracked, ', newline num2str(dN_cutoff(1)*100), ...
            '$\;<d_{\rm{NU}}<\;$', num2str(dN_cutoff(2)*100), 'cm']})
end

nexttile
hold on; grid on; grid minor
plot(IPop(:,2),IPop(:,ind_B),'-k*')
plot(IPop_interp(:,1),IPop_interp(:,ind_B),'-r',LineWidth=1.2)
% ylabel('Rocket bodies [obj]'); xlabel('Altitude [km]'); 
title('Rocket bodies')
xlim([IPop(1,2) IPop(end,2)])

% han=axes(fig,'visible','off');
% han.Title.Visible='on'; 
% title(han,['Initial population, h$_{shell} =$ ', num2str(h_bin)])
leg = legend('MASTER','Interpolation','Orientation', 'Horizontal');
leg.Layout.Tile = 'north';
title(tl, ['Initial population, h$_{shell} =$ ', num2str(h_bin)],'interpreter','latex','fontsize', 24);
ylabel(tl,'Number of objects [-]','interpreter','latex','fontsize', 22);
xlabel(tl,'Altitude [km]','interpreter','latex','fontsize', 22);




