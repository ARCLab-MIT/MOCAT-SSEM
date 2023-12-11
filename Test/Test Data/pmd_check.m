% pmd_check

function my_sim = pmd_check(pmd, disposal_alt)
% This function takes a mass, radius, starting altitude, and a
% goal_lifetime and finds the closest altitude bun that will deorbit in
% less than the goal_lifetime.

%% Make Scenario Properties

%goal_lifetime = 25

fileName = matlab.desktop.editor.getActiveFilename;
folder = fileparts(which(fileName)); 
addpath(genpath(folder));
cd(folder);
addpath(genpath(pwd))
addpath(genpath(folder));

disp("Pm: " + string(pmd) + " Disposal altitude" + string(disposal_alt));

%% Scenario Properties

% Add the time stuff needed
VAR = MOCATSSEM_VAR_Cons();
VAR.start_year = 2022;
VAR.start_month = 12;
VAR.start_day = 1;
VAR.start_hour = 0;
VAR.start_minute = 0;
VAR.start_second = 0;
VAR.start_second_fraction = 0;
VAR.simulation_duration = 200.0; %Years
VAR.N_steps = 200;
VAR.scen_times = linspace(0,VAR.simulation_duration,VAR.N_steps)';
VAR.density_filepath = fullfile(".", "Atmosphere Model", "JB2008", "Precomputed", "dens_highvar_2000.mat");
VAR.dens_model = @JB2008_dens_func;
VAR.time_dep_density = true;

%VAR.dens_model = @static_exp_dens_func;
%VAR.time_dep_density = false;

VAR.integrator = @ode15s; % Fixed options: ode1, ode2, ode3, ode4, ode5
%                                     % Variable options: ode45, ode21,
%                                     % ode113, oder78, ode89, ode15s,
%                                     % oder23s, ode23t, ode23tb
VAR.integrator_type = 'variable'; %or 'fixed';

scen_properties = scen_properties_class(VAR);



%% Make Species

% S
% species_properties = {};
% species_properties.sym_name = "S";
% species_properties.Cd = 2.2; % unitless
% species_properties.mass = 1250; % kg
% species_properties.radius = 4.1; % m
% species_properties.A = 24; % m^2 %TODO don't guess area.
% species_properties.amr = species_properties.A./species_properties.mass;  % m^2/kg
% species_properties.beta = species_properties.Cd*species_properties.amr; % ballistic coefficient
% species_properties.drag_effected = false; % bool, does object decrease altitude due to drag
% species_properties.active = true; % bool, use alpha_active for collisions vs. trackable
% species_properties.maneuverable = true; % bool, use alpha_active for collisions vs. trackable
% species_properties.trackable = true; % bool, others can avoid this with alpha.
% species_properties.deltat = 5.0; % lifetime of spacecraft in years
% species_properties.Pm = pmd; % Post-mission disposal efficacy Float [0,1]
% species_properties.alpha = 1e-4; % Efficacy of Collision Avoidance vs. inactive
% species_properties.alpha_active = 1e-4; % Efficacy of Collision Avoidance vs. other active 
% species_properties.RBflag = 0; % No Rocket Bodies yet.
% species_properties.lambda_constant = 2000;
% S_species = species(@launch_func_constant, @pmd_func_sat, @drag_func, species_properties, scen_properties);

% S
species_properties.sym_name = "S";

% Checked properties
species_properties.Pm = pmd; % Post-mission disposal efficacy Float [0,1]
species_properties.disposal_altitude = disposal_alt;

% Rest of properties.
species_properties.Cd = 2.2; % unitless
species_properties.mass = [255]; % kg
species_properties.radius = [sqrt(15.45/pi)]; % m
species_properties.A = species_properties.radius.^2 * pi; % m^2 %TODO don't guess area.
% For this amr and others, we use first mass so multi_property_species will
% scale it automatically for us using the scalred radii.
species_properties.amr = species_properties.A./species_properties.mass;  % m^2/kg
species_properties.beta = species_properties.Cd*species_properties.amr; % ballistic coefficient

% Capabilities
species_properties.drag_effected = false; % bool, does object decrease altitude due to drag
species_properties.active = true; % bool, use alpha_active for collisions vs. trackable
species_properties.maneuverable = true; % bool, use alpha_active for collisions vs. trackable
species_properties.trackable = true; % bool, others can avoid this with alpha.
species_properties.deltat = [5]; % lifetime of spacecraft in years
species_properties.alpha = 1e-5; % Efficacy of Collision Avoidance vs. inactive
species_properties.alpha_active = 1e-5; % Efficacy of Collision Avoidance vs. other active 
species_properties.RBflag = 0; % No Rocket Bodies yet.
species_properties.lambda_constant = 2000;
S_species = species(@launch_func_constant, @pmd_func_sat, @drag_func, species_properties, scen_properties);
S_masses = S_species.species_properties.mass;
S_radii = S_species.species_properties.radius;

% N
species_properties.sym_name = "N";
species_properties.Cd = 2.2; % unitless
species_properties.deltat = NaN;
species_properties.disposal_altitude = NaN;
species_properties.mass =   [0.001418952683 1.418952683 15  S_masses]; 
species_properties.radius = [0.005          0.05        0.1097348193 S_radii]; % m
species_properties.A = pi*species_properties.radius.^2;
species_properties.amr = species_properties.A./species_properties.mass;  % m^2/kg
species_properties.beta = species_properties.Cd*species_properties.amr; % ballistic coefficient
species_properties.drag_effected = true; % bool, does object decrease altitude due to drag
species_properties.maneuverable = false; % bool, use alpha_active for collisions vs. trackable
species_properties.active = false;
N_species = multi_property_species( @launch_func_null, @pmd_func_derelict,@drag_func, species_properties, scen_properties).species_list;

N_species = pair_actives_to_debris(scen_properties,S_species, N_species); % Also adjusts handles, so re-assignment is not strictly necessary

scen_properties.species = [S_species N_species];
scen_properties.species_cell = struct('S', S_species, 'N', N_species);

%[scen_properties] = make_collision_pairs_SBM(scen_properties);

% make x0
alt_bin = find_alt_bin(disposal_alt, scen_properties);
x0 = zeros(scen_properties.N_shell * numel(scen_properties.species), 1);

%%
my_sim = simulation_class(scen_properties.species, scen_properties);
my_sim.build_model();
my_sim.scen_properties.integrator = @ode15s;
my_sim.scen_properties.disp_times = false;

%%
my_sim.run_model(x0, 'progressBar', true, 'disp_times', false);
my_sim.pop_vs_time_vis()
my_sim.total_species_evol_vis()
disp("Done")
