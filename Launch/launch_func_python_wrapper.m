function [Lambdadot] = launch_func_python_wrapper(x, t, h, h_i, species_properties, scen_properties)
    %   Wrapper that tries to improve efficency by limiting the number of
    %   python calls. This will ensure that at each state X, and time T,
    %   the python call only occurs once.  The code will then return the
    %   proper lambda value from the cached Lambdadot variable based on the
    %   relevant shell.
    %   x is the state of the system during integration
    %   t is time from scenario start in years
    %   h is the set of altitudes of the scenario above ellipsoid in km of shell lower edges (unused).
    %   h_i is the index of the particular altitude shell in question.
    %   species_properties is a structure with properties for the species
    %   scen_properties is a structure with properties for the scenario
    %   Lambdadot is the rate of change in the species (sats/year) in the given shell at the specified time due to launch.
    %TODO: check hmid

    % If the X state changes or there is a time update, call the python method
    % again.
    if any(x ~= species_properties.last_calc_x )|| any(t ~= species_properties.last_calc_t)
        
        % Try to do a workaround saving and loading mat files
        % species_properties.lambda_python_args.name = species_properties.sym_name + "_launch";
        % name =  py.importlib.import_module('test_launch_fun').python_launch_save(x, t, h, species_properties.lambda_python_args); %py.test_launch_fun.python_launch_save(x, t, h, species_properties.lambda_python_args).';
        % load(species_properties.lambda_python_args.name, "Lambdadot");
        
        % Type cast from python to matlab
        if ~exist('test_launch_fun', 'var')
            test_launch_fun = py.importlib.import_module('test_launch_fun');
        end
        
        Lambdadot = test_launch_fun.python_launch_func2(x, t, h, species_properties.lambda_python_args).';
        Lambdadot = double(Lambdadot).';
        
        species_properties.last_calc_x = x;
        species_properties.last_calc_t = t;
        species_properties.last_lambda_dot = Lambdadot;
    end

    % Return lambda at the specific altitude of interest using the existing or updated launch function.
    Lambdadot = species_properties.last_lambda_dot(h_i);

    end
