function [Lambdadot] = launch_func_model_loop_builder(t, h, species_properties, scen_properties)
    % Function to use to insert placeholders in full_lambda during build.
    %   t is time from scenario start in years
    %   h is the set of altitudes of the scenario above ellipsoid in km of shell lower edges.
    %   species_properties is a structure with properties for the species
    %   scen_properties.launch_filepath = relative or absolute path to
    %   launch csv file with column for time (year) and one or more launch
    %   rates
    %   scen_properties is a structure with properties for the scenario
    %   Lambdadot is the rate of change in the species in each shell at the specified time due to launch.

    % Add this to the list of functions that get called every loop of the
    % ODE operator. Don't need to do anything with the output because it
    % modifies scen_properties.Lambdadot directly.
    calling_func = dbstack(1).name;
    if  calling_func == "simulation_class.build_model"
        func_name = str2func("launch_func_model_loop");
        launch_func_model_loop_placeholder_func = @(x, t) 0;
        scen_properties.functions_to_run_each_model_loop{end+1} = @(x,t) func_name(t, h, species_properties, scen_properties);
        Lambdadot = repelem({@(x,t) launch_func_model_loop_placeholder_func(x,t)}, scen_properties.N_shell).';
    end