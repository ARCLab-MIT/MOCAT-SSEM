classdef sim_optimizer < handle
    properties
        simulation %instance of simulation class for the relevant optimization
        simulation_run_settings
        
        objective %function with signature objective(simulation) that takes 
        % the results from the run simulation, calculates the objective
        % value, and returns it.
        objective_settings
        
        optimizer
        optimizer_settings
        input_adjust_sim

        x0
        optvar_min % Min value for each optimization variable
        optvar_max % Max value for each optimization variable
        N_var % Number of optimization vars, set by 

        % Generated
        opj_func % Takes build optimization model

        % Results
        optimizer_results

    end % end properties
    methods
        function obj = sim_optimizer(simulation, objective, optimizer, ...
                                     input_adjust_sim, optimizer_settings, ...
                                     objective_settings, ...
                                     simulation_run_settings, ...
                                     x0, optvar_min, optvar_max, N_var)
            obj.simulation = simulation; %simulation_class object
            obj.objective = objective;
            obj.optimizer = optimizer;
            obj.input_adjust_sim = input_adjust_sim;
            %obj.sim_wrapper_function = sim_wrapper_function;

            % Optional
            % TODO Check if parameters are passed.
            obj.optimizer_settings = optimizer_settings;
            obj.objective_settings = objective_settings;
            obj.simulation_run_settings = simulation_run_settings;
            obj.x0 = x0;
            obj.optvar_min = optvar_min;
            obj.optvar_max = optvar_max;
            obj.N_var = N_var;

            % wrapper function

            obj.opj_func = @(lambda)wrap_sim(obj, lambda);

        end

        function [J] = wrap_sim(obj, lambda)
            
            obj.input_adjust_sim(obj.simulation, lambda);
            obj.simulation.run_model(obj.x0, obj.simulation_run_settings) 
            [J, results] = obj.objective(obj.simulation, obj.objective_settings);
            %disp("J = " + num2str(J, '%.e '))
        end

        function [xg_monte,Jg_monte,exitflag,output] = run_optimization(obj)
            tic
            [xg_monte,Jg_monte,exitflag,output] = obj.optimizer( ...
                                                                obj.opj_func, ...
                                                                obj.N_var, ...
                                                                obj.optvar_min, ...
                                                                obj.optvar_max, ...
                                                                obj.optimizer_settings);       
            toc
            obj.optimizer_results = {xg_monte,Jg_monte,exitflag,output};
        end
    end %end methods
end %end class