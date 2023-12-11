classdef species < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        launch_func
        pmd_func
        drag_func
        species_properties
        scen_properties
    end

    methods
        function obj = species(launch_func, pmd_func, drag_func, species_properties, scen_properties)
            %UNTITLED Construct an instance of this class
            % launch_func: a function that specifies launch rate. Launch 
            %   functions takes a parameter t (float), a species-level 
            %   properties (structure), and a model -level properties 
            %   structure. It returns the launch quantity at the time for 
            %   this species as a function of time.
            %   TODO: figure out of this will be symbolic of quantity.
            % pmd_func: a function that specifies the PMD rate. The PMD 
            %   functions takes a parameter t (float), a species-level 
            %   properties (structure), and a model -level properties 
            %   structure. It returns the PMD rate for this species as a
            %   function of time.
            %   TODO: figure out of this will be symbolic of quantity.
            % drag_func: a function that specifies the reduction in 
            %   altitude per unit time. The drag functions takes a 
            %   parameter t (float), h (altitude in km), a species-level 
            %   properties (structure), and a model -level properties 
            %   structure. It returns the PMD rate for this species as a
            %   function of time.
            obj.launch_func = launch_func;
            obj.pmd_func = pmd_func;
            obj.drag_func = drag_func;
            species_properties.sym_name = strrep(species_properties.sym_name, ".", "p"); %avoid . in symbolic var name
            species_properties.sym = sym(species_properties.sym_name + "_",[scen_properties.N_shell,1]);
            obj.species_properties = species_properties_class(species_properties);
            obj.scen_properties = scen_properties;
        end
    end
end

% launch_func(t, species_properties, scen_properties)
% pmd_func(t, species_properties, scen_properties)
% drag_func(t, h, species_properties, scen_properties)

% assemble model.
% add symbolic variable to simulation.

% species_properties = {};
% species_properties.Cd = 2.2; % unitless
% species_properties.mass = 223.0; % kg
% species_properties.radius = 0.5; % m
% species_properties.A = 1.741; % m^2
% species_properties.amr = species_properties.A/species_properties.radius;  % kg
% species_properties.slotted = false; % bool
% species_properties.slotting_effectiveness  = 1.0; % Float [0,1],
% species_properties.maneuveable = true; % bool, use alpha_active for collisions vs. trackable
% species_properties.trackable = true; % bool, others can avoid this with alpha.
% species_properties.deltat = 5.0; % lifetime of spacecraft in years
% species_properties.Pm = .95; % Post-mission disposal efficacy Float [0,1]
% species_properties.alpha = 0.2; % Efficacy of Collision Avoidance vs. inactive
% species_properties.alpha_active = 0.001; % Efficacy of Collision Avoidance vs. other active 
% species_properties.orbit_raising = 'onstation'; % Whether orbit-raising will occur.
% species_properties.insert_altitude = 250.0; % Altitude in km for insertion altitude
% species_properties.onstation_altitude = 7500.0; % Altitude in km for onstation altitude
% 

