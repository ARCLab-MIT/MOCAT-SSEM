function [Lambdadot] = launch_func_model_loop(t, h, species_properties, scen_properties)
    % Updates scen_properties.full_lambda and returns Lambdadot for the
    % chosen species. This is the function that gets run at each ODE
    % iteration.
    %   t is time from scenario start in years
    %   h is the set of altitudes of the scenario above ellipsoid in km of shell lower edges.
    %   species_properties is a structure with properties for the species
    %   scen_properties.launch_filepath = relative or absolute path to
    %   launch csv file with column for time (year) and one or more launch
    %   rates
    %   scen_properties is a structure with properties for the scenario
    %   Lambdadot is the rate of change in the species in each shell at the specified time due to launch.
    
    %disp("Main_Scenario: " + string(scen_properties.isMain == true) + ", launch_func_model_loop called for species " + ...
    %     species_properties.sym_name + " at time " + num2str(t) + ...
    %     " with t_plan_max " + num2str(species_properties.t_plan_max) + ".") 

    if  t > species_properties.t_plan_max && scen_properties.isMain == true
        % Load precompiled model
        load(scen_properties.saved_model_path, ["my_sim"]);
        my_sim.scen_properties.isMain = false;

        % Run model to anticipate future step, update planning horizon
        species_properties.t_plan_max = t + species_properties.t_plan_period;
        disp("Updating t_plan_max to " + num2str(species_properties.t_plan_max) + ".") 
        start_date = scen_properties.start_date + years(t);
        disp("Running secondary simulation starting for species " + ...
            species_properties.sym_name + " from " + string(start_date) + ...
            " and continuing for " + ...
            num2str(species_properties.t_plan_period) + " years.") 
        my_sim.run_model(scen_properties.X, 'progressBar', false, 'start_date', start_date, 'simulation_duration',species_properties.t_plan_period, 'N_step', 5);
        disp("Completed secondary simulation.") 
        species_properties.prev_prop_results = my_sim.results.X(end,my_sim.scen_properties.species(1).species_properties.eq_idxs).';
        
        % Modify launch rate based on it
        Lambdadot = cell(size(species_properties.eq_idxs.'));
        launch_func_model_lambdas = species_properties.prev_prop_results + 10*ones(size(species_properties.eq_idxs.'));
        for ii = 1:length(launch_func_model_lambdas)
            Lambdadot{ii} = @(x, t) launch_func_model_lambdas(ii, 1);
        end 
        scen_properties.full_lambda(species_properties.eq_idxs, 1) = Lambdadot;
        disp("New Lambda for " + species_properties.sym_name + ":")
        disp(cellfun(@(fun) fun(species_properties.prev_prop_results, t),Lambdadot,'UniformOutput',false).');
        disp("")
    else %Being called by ODE builder, but within plan period.
        Lambdadot = scen_properties.full_lambda{species_properties.eq_idxs, 1};
    end

end
