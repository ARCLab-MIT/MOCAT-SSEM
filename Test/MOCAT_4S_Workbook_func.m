function [my_sim] = MOCAT_4S_Workbook_func()
% Function with the same code as MOCAT 4S Workbook

% The demo workbooks are in a subfolder of the main directory. We need to
% add the path recursively for the main directory.
fileName = matlab.desktop.editor.getActiveFilename;
[folder, ~, ~] = fileparts(fileName);
folder = fullfile(folder, '..');
addpath(fileparts(folder));
addpath(genpath(folder));
cd(folder);
scenario_properties = MOCATSSEM_Scen_Prop_Cons( ...
    start_date = datetime(2022,1,30,9,42,28), ... # Starting date (year, month, day, hour, minute, second)
    simulation_duration = 100.0, ... % Years of simulation to run
    steps = 1000, ... % Number of steps to run in simulation
    min_altitude = 200, ... % Maximum altitude shell in km
    max_altitude = 900, ... % Minimum altitude shell in km
    shells = 5, ... % Number of altitude shells
    ... % We create this many equally spaced shells between the minimum and maximum altitudes
    delta = 10, ... %ratio of the density of disabling to lethal debris
    ... % this term, delta, considers the possibility that disabling collisions can generate new derelicts
    integrator = @ode15s); % MATLAB integrator used to run simulation
    % we typically use ode15s because often, models with many shells become stiff diff eqs.


atm_density_choice = "Static Exponential";
switch atm_density_choice
    case "JB2008"
        scenario_properties.density_filepath = fullfile(".","Drag","Atmosphere Model", "JB2008", "Precomputed", "dens_highvar_2000.mat");
        scenario_properties.dens_model = @JB2008_dens_func;
        scenario_properties.time_dep_density = true;
    case "Static Exponential"
        scenario_properties.dens_model = @static_exp_dens_func;
        scenario_properties.time_dep_density = false;
end
scen_properties = scen_properties_class(scenario_properties);
scen_properties.options.AbsTol = 1e-2;
scen_properties.options.RelTol = 1e-2;
% scen_properties.options
% Define Species Properties
launch_func_choice = "Null";
switch launch_func_choice
    case "Null"
        launch_func_sat = @launch_func_null;
    case "Constant"
        launch_func_sat = @launch_func_constant;
        lambda_constant = num2cell(500 * rand(scen_properties.N_shell,1));
end
my_drag_func = @drag_func;
species_properties = struct;

alpha_active = 0.01; % Efficacy of Collision Avoidance vs. other active 
alpha = 2e-3; % Efficacy of Collision Avoidance vs. inactive

% Su
species_properties.sym_name = "Su";
species_properties.Cd = 2.2; % unitless
species_properties.mass = 223; % kg
species_properties.radius = 1.490/2; % m
species_properties.A = 1.741; % m^2
species_properties.amr = species_properties.A/species_properties.mass;  % m^2/kg
species_properties.amr = species_properties.A./species_properties.mass;  % m^2/kg
species_properties.beta = species_properties.Cd*species_properties.amr; % ballistic coefficient

% Orbit Properties
species_properties.slotted = false; % bool

% Capabilities
species_properties.drag_effected = false; % bool, does object decrease altitude due to drag
species_properties.active = true; % bool, use alpha_active for collisions vs. trackable
species_properties.maneuverable = true; % bool, use alpha_active for collisions vs. trackable
species_properties.trackable = true; % bool, others can avoid this with alpha.
species_properties.deltat = 5.0; % lifetime of spacecraft in years
species_properties.Pm = .90; % Post-mission disposal efficacy Float [0,1]
species_properties.alpha = 2e-3; % Efficacy of Collision Avoidance vs. inactive
species_properties.alpha_active = 1e-3; % Efficacy of Collision Avoidance vs. other active 
species_properties.RBflag = 0; % No Rocket Bodies yet.
if launch_func_choice == "Constant"
    species_properties.lambda_constant = lambda_constant; 
end
Su_species = species(launch_func_sat, @pmd_func_sat, my_drag_func, species_properties, scen_properties);
species_properties = struct;

% S
species_properties.sym_name = "S";
species_properties.Cd = 2.2; % unitless
species_properties.mass = 223; % kg
species_properties.radius = 1.490/2; % m
species_properties.A = 1.741; % m^2
species_properties.amr = species_properties.A/species_properties.mass;  % m^2/kg
species_properties.amr = species_properties.A./species_properties.mass;  % m^2/kg
species_properties.beta = species_properties.Cd*species_properties.amr; % ballistic coefficient

% Orbit Properties
species_properties.slotted = true; % bool
species_properties.slotting_effectiveness  = 1.0; % Float [0,1],

% Capabilities
species_properties.drag_effected = false; % bool, does object decrease altitude due to drag
species_properties.active = true; % bool, use alpha_active for collisions vs. trackable
species_properties.maneuverable = true; % bool, use alpha_active for collisions vs. trackable
species_properties.trackable = true; % bool, others can avoid this with alpha.
species_properties.deltat = 5.0; % lifetime of spacecraft in years
species_properties.Pm = .90; % Post-mission disposal efficacy Float [0,1]
species_properties.alpha = 2e-3; % Efficacy of Collision Avoidance vs. inactive
species_properties.alpha_active = 1e-3; % Efficacy of Collision Avoidance vs. other active 
species_properties.RBflag = 0; % No Rocket Bodies yet.
if launch_func_choice == "Constant"
    species_properties.lambda_constant = lambda_constant; 
end

S_species = species(launch_func_sat, @pmd_func_sat, my_drag_func, species_properties, scen_properties);
species_properties = struct;

% D
species_properties.sym_name = "D";
species_properties.Cd = Su_species.species_properties.Cd;
species_properties.mass = Su_species.species_properties.mass;
species_properties.radius = Su_species.species_properties.radius;
species_properties.A = Su_species.species_properties.A;
species_properties.amr = Su_species.species_properties.A/Su_species.species_properties.mass(1);  % m^2/kg
species_properties.beta = Su_species.species_properties.Cd*Su_species.species_properties.amr; % ballistic coefficient
species_properties.drag_effected = true; % bool, does object decrease altitude due to drag
species_properties.maneuverable = false; % bool, use alpha_active for collisions vs. trackable
species_properties.active = false;

D_species = species(@launch_func_null, @pmd_func_derelict, my_drag_func, species_properties, scen_properties);
D_species.species_properties.pmd_linked_species = [S_species, Su_species];
species_properties = struct;

% N
species_properties.sym_name = "N";
species_properties.Cd = 2.2; % unitless
species_properties.mass = 0.640; 
species_properties.radius = 0.180/2; % m
species_properties.A = 0.020;
species_properties.amr = species_properties.A/species_properties.mass;  % m^2/kg
species_properties.beta = species_properties.Cd*species_properties.amr; % ballistic coefficient
species_properties.drag_effected = true; % bool, does object decrease altitude due to drag
species_properties.maneuverable = false; % bool, use alpha_active for collisions vs. trackable
species_properties.active = false;

N_species = species(@launch_func_null, @pmd_func_none, my_drag_func, species_properties, scen_properties);
species_properties = struct;

scen_properties.species = [S_species, D_species, N_species, Su_species];
scen_properties.species_cell = struct('S', S_species,'D', D_species, 'N', N_species, 'Su', Su_species);



% Collision Modeling
% The make_collision_pairs_MOCAT4.m function calculates collision coefficients for the table above, based on the species properties and NASA's standard break-up model, and adds the collision terms to the source sink equations. The make_collision_pairs_MOCAT4.m function is specific to MOCAT-4. The generalized function for collision modelling for other MOCAT versions is make_collision_pairs_SBM.m
scen_properties = make_collision_pairs_MOCAT4(scen_properties, Su_species, S_species, D_species, N_species);

% Define Initial Population
% This section creates an initial population for satellites, derelict and debris. The MC2SSEM_population function takes a list of space objects from initialized.mat and categorizes them into satellites, derelict and debris. x0 is a total count of each species in each shell that we will start the simulation from.
% Species order and definition (NOTE: do not change the order): 

scen_properties.species_types = [isfield(scen_properties.species_cell,'S') && ~isempty(scen_properties.species_cell.S),...
                                 isfield(scen_properties.species_cell,'D') && ~isempty(scen_properties.species_cell.D),...
                                 isfield(scen_properties.species_cell,'N') && ~isempty(scen_properties.species_cell.N),...
                                 isfield(scen_properties.species_cell,'Su') && ~isempty(scen_properties.species_cell.Su),...
                                 isfield(scen_properties.species_cell,'B') && ~isempty(scen_properties.species_cell.B),...
                                 isfield(scen_properties.species_cell,'U') && ~isempty(scen_properties.species_cell.U)];
load('initialized.mat', 'sats');
x0 = MC2SSEM_population(sats,scen_properties);

% Display starting population
for spec_index = 1:length(scen_properties.species)
   spec_names(spec_index) = scen_properties.species(spec_index).species_properties.sym_name;
end

T = array2table(x0, 'VariableNames', spec_names);
disp(array2table([[1:scen_properties.N_shell]',x0], 'VariableNames', ['Shell', spec_names]))

% Build Model
% simulation_class is responsible for building and running the model. The build_model() function compiles the equations. It adds together all of the terms (regarding launch, collisions, post-mission disposal, and atmospheric drag).
% Building the model is only required once in order to run the model if the user doesn't want to change any input parameters.
x0 = reshape(table2array(T), [], 1);% Pad initial conditions with the indicator variables
my_sim = simulation_class(scen_properties.species, scen_properties);
my_sim.build_model();

% Run Simulation
% run_model() is where the simulation model is actually run. It takes a list of initial conditions (x0).
% The function also accepts many optional parameters: 
x0 = reshape(x0, [], 1);
my_sim.run_model(x0, 'progressBar', true, 'disp_times', false) 

% % Plots
% % Produces a Figure with subplots for each species population evolution over time 
% my_sim.total_species_evol_vis();
% % Produces a Figure with subplots for each species population density evolution over time 
% my_sim.density_species_evol_vis();
% % Produces a Figure for the evolution of the total population of each species over time across all shells 
% my_sim.pop_vs_time_vis();
% % Graphs the rate of change for the total inactive population added across all altitudes at each time. 
% % percentage = (bool, defaults to false) that expresses rate of change in percentage terms. 
% my_sim.total_deb_pop_time_deriv('percentage', true) 
% % Produces a Figure with subplots for numerical derivative of each inactive species population evolution over time constraint is the value used in the color scale.
% % percentage is bool (default false) for whether to use percentage or absolute.
% % color_scale_crop_percentiles is pair (default [0, 100]). Color scale will be  set based on the inclusive range from the lower value to the higher value.
% my_sim.total_deb_species_deriv_evolv_vis();
% % Produces a Figure for the evolution of the total population of all objects over time across all shells 
% my_sim.tot_vs_time_vis();
% % Produces a Figure per species for the evolution of the total population of each species over time across all shells 
% my_sim.species_time_vis();
% % Displays total spatial density in each shell over time 
% my_sim.density_vs_time();
% % Computes and displays Number_Time_Product
% my_sim.num_time_prod(x0);
% % Produces a Figure comparing the total initial number of objects and the total final number of objects 
% my_sim.init_vs_final_vis();
% % Produces a Figure comparing the initial population of each species alongside a plot of the final population of each species 
% my_sim.init_and_final_vis();
% % Produces a Figure showing number of debris across altitudes by debris type 
% my_sim.final_debris_vis();
% % Computes and displays cumulative Criticality of Space Index - CSI
% my_sim.cum_CSI()
% Close all figures created (or run this line in your command line)
% close all

end
