function [VAR] = MOCATSSEM_VAR_Cons()

warning("This method is being deprecated. Please switch your workflow to use MOCATSSEM_Scen_Prop_Cons") 

% Constants
rad = pi/180;
mu_km=398600.4415;%km^3/s^2
mu=3.986004418e14;%meters^3/s^2
VAR.re = 6378.1366; % [km]
years_s = 365*24*3600;

% Parameters needed for all MOCAT-SSEMs
VAR.N_shell = 40;
VAR.h_max = 1400;
VAR.h_min = 200;
VAR.N_step = 5; % how many steps to calculate between t0 and tf
R0=linspace(VAR.h_min,VAR.h_max,VAR.N_shell+1);
R02=linspace(VAR.h_min,VAR.h_max,VAR.N_shell+1);
VAR.HMid = R02(1:end-1) + diff(R02)/2;
VAR.deltaH=R02(2)-R02(1); % thickness of the shell [km]

%VAR.v_imp = 10; % impact velocity [km/s]
R1=R0;
R0=(VAR.re+R0)*1000;
% V=4*pi*R0.^2*deltaH*1000;
VAR.V=4/3*pi*(diff(R0.^3)); % volume of the shells [m^3]
VAR.v_imp2 = 10 * ones(size(VAR.V)); % impact velocity [km/s] %Shell-wise
VAR.LC = 0.01; % minimum size of fragments [m]

VAR.v=VAR.v_imp2*1000*(24*3600*365.25);% impact velocity [m/year]
%VAR.v2=VAR.v_imp*1000*(24*3600*365.25);% impact velocity [m/year]
VAR.Dhl=VAR.deltaH*1000;
VAR.Dhu=-VAR.deltaH*1000;

%Integration stuff
options = odeset('reltol', 1.e-4,'abstol', 1.e-4);
VAR.options = options;
VAR.mu = mu;
VAR.R0 = R0;
VAR.R02 = R02;

%launch_rate = 'generic';
%VAR.S_Su = 0.1; %ratio of unslotted to slotted spacecraft with generic or distribution model
%MOCAT = '4S'; 
%VAR.launch_rate = launch_rate;



