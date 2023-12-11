function [Lambdadot] = launch_func_python(t, h, species_properties, scen_properties)
    % Adds constant launch rate from species_properties.lambda_constant
    %   t is time from scenario start in years
    %   h is the set of altitudes of the scenario above ellipsoid in km of shell lower edges.
    %   species_properties is a structure with properties for the species
    %   scen_properties.launch_filepath = relative or absolute path to
    %   launch csv file with column for time (year) and one or more launch
    %   rates
    %   scen_properties is a structure with properties for the scenario
    %   Lambdadot is the rate of change in the species in each shell at the specified time due to launch.
    
    if ~isfield(species_properties, "lambda_python_args")
        species_properties.lambda_python_args = {};
    end

    tf = isfield(species_properties.lambda_python_args, "idx");
    if ~tf
        [~,idx] = intersect(string(species_properties.sym),string(scen_properties.var_col),'stable');
        species_properties.lambda_python_args.idx = idx;
    end

    %x = scen_properties.X;
    lambda_python_args = species_properties.lambda_python_args;
    
    Lambdadot = cell(size(h,2), 1);
    for h_i = 1:size(h,2)
        lambda_python_args.cur_idx = species_properties.lambda_python_args.idx(h_i);
        lambda_python_args.col_idx = scen_properties.indicator_var_list(2).indicator_idxs(lambda_python_args.cur_idx);
        Lambdadot{h_i, 1} = @(x, t) double(py.test_launch_fun.python_launch_func(x, t, h(h_i), lambda_python_args));
    end
    end