% drag_check

function my_sim = applicant_drag_check(scen_properties, mass, radius, starting_altitude, candidate_cycle_launches)
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

%% Species Properties
species_properties = {};

% Capabilities
species_properties.drag_effected = false; % bool, does object decrease altitude due to drag
species_properties.active = true; % bool, use alpha_active for collisions vs. trackable
species_properties.maneuverable = true; % bool, use alpha_active for collisions vs. trackable
species_properties.trackable = true; % bool, others can avoid this with alpha
species_properties.RBflag = 0; % No Rocket Bodies yet.


%% Make Species
t = 10;
my_drag_func = @drag_func;
null_launch_func = @launch_func_null;
launch_func = @launch_func_gauss;
pmd_func = 0;

% C
% Assuming al spheres of 1 cm, 10cm diameter, 15 kg 
species_properties.sym_name = "C";
species_properties.Cd = 2.2; % unitless
species_properties.deltat = NaN;
species_properties.mass =   [mass]; 
species_properties.radius = [radius]; % m
species_properties.A = pi*species_properties.radius.^2;
species_properties.amr = (species_properties.A)./species_properties.mass;  % m^2/kg
species_properties.beta = species_properties.Cd*species_properties.amr; % ballistic coefficient
species_properties.drag_effected = true; % bool, does object decrease altitude due to drag
species_properties.maneuverable = false; % bool, use alpha_active for collisions vs. trackable
species_properties.active = false;
species_properties.lambda_constant = zeros(scen_properties.N_shell, 1);
species_properties.lambda_constant(find_alt_bin(starting_altitude, scen_properties)) = species_properties.lambda_constant(find_alt_bin(starting_altitude, scen_properties));
C_species = species(@launch_func_null, @pmd_func_none,my_drag_func, species_properties, scen_properties);


scen_properties.species = C_species;
scen_properties.species_cell = struct('C', C_species);

file_path = candidate_cycle_launches;
[x0,FLM_steps] = ADEPT_Traffic_Model(file_path, scen_properties);
FLM_to_launch_functions(FLM_steps, scen_properties);

% make x0
alt_bin = find_alt_bin(starting_altitude, scen_properties);
x0 = zeros(scen_properties.N_shell, 1);
x0(alt_bin) = 1;
%%
my_sim = simulation_class(scen_properties.species, scen_properties);
my_sim.build_model();
my_sim.scen_properties.integrator = @ode15s;
my_sim.scen_properties.disp_times = false;


%%
my_sim.run_model(x0, 'progressBar', true, 'disp_times', false);
%my_sim.pop_vs_time_vis()
%my_sim.total_species_evol_vis()
%disp("Done")
