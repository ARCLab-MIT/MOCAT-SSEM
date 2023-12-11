classdef scen_properties_class < handle
    properties
        re
        N_shell
        h_max
        h_min
        N_step
        deltaH
        alpha
        alpha_active
        Dt
        slot_par
        delta
        P
        Cd
        mass
        A
        diameter
        v_imp
        v_imp2
        LC
        deg_intrinsic
        sep_dist_method
        sep_angle
        sep_dist
        input_pop
        V
        v
        Dhl
        Dhu
        area_mass
        x0
        options
        beta
        mu
        K0
        phi
        R0
        R02
        HMid
        N_sat
        S_Su
        launch_rate
        
        % Set for scenario
        % Density information
        density_filepath
        dens_model
        interp_dens_var % Set during build_model(obj) based on times
        time_dep_density
        
        sym_drag
        full_drag
        sym_lambda
        full_lambda

        start_year
        start_month
        start_day
        start_hour
        start_minute
        start_second
        start_second_fraction
        start_date
        simulation_duration
        N_steps
        scen_times

        species
        species_pairs
        % Define the species considered and the initial population
        % species:[S,D,N,Su,B,U];
        % S = Active slotted satellites
        % D = Derelicts
        % N = Debris objects
        % Su = Active unslotted satellites
        % B = Rocket bodies
        % U = Lethal non-trackable debris objects
        species_types
        species_cell % Struct with S, D, N, Su, B arrays or whatever species types exist

        integrator
        integrator_type

        integrated_indicator_var_list
        num_integrated_indicator_vars = 0
        indicator_idxs 

        indicator_var_list

        X % Used to understand current state of the system for python lauunch functions
        t % Used to understand current state of the system for python lauunch functions
        var_col

        saved_model_path
        functions_to_run_each_model_loop = {}
        isMain = true % is the main simulation vs. a nested side simulation.

        disp_times %bool for whether to print time steps during integration

    end % end properties
    methods
        function obj = scen_properties_class(struct)
            for field = fieldnames(struct)'    %enumerat fields
                try
                  obj.(field{1}) = struct.(field{1});   %and copy
                catch
                  warning('Error copying field %s', field{1});
                end
            end
            obj.start_date = datetime(obj.start_year, ...
                                      obj.start_month, ...
                                      obj.start_day, ...
                                      obj.start_hour, ...
                                      obj.start_minute, ...
                                      obj.start_second, ...
                                      obj.start_second_fraction);

        end
    end %end methods
end %end class
