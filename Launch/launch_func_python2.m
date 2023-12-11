function [Lambdadot] = launch_func_python2(t, h, species_properties, scen_properties)
    % Adds constant launch rate from species_properties.lambda_constant
    %   t is time from scenario start in years
    %   h is the set of altitudes of the scenario above ellipsoid in km of shell lower edges.
    %   species_properties is a structure with properties for the species
    %   scen_properties.launch_filepath = relative or absolute path to
    %   launch csv file with column for time (year) and one or more launch
    %   rates
    %   scen_properties is a structure with properties for the scenario
    %   Lambdadot is the rate of change in the species in each shell at the specified time due to launch.
    
    if ~isprop(species_properties, "lambda_python_args")
        species_properties.lambda_python_args = {};
    end

    if ~isfield(species_properties.lambda_python_args, "idx")
        species_properties.lambda_python_args.idx = uint16(species_properties.eq_idxs);
    end

    if ~isfield(species_properties.lambda_python_args, "col_idx")
        species_properties.lambda_python_args.col_idx = uint16(scen_properties.indicator_var_list(2).indicator_idxs);
    end

    Lambdadot = cell(size(h,2), 1);
    for h_i = 1:size(h,2)
        Lambdadot{h_i, 1} = @(x, t) launch_func_python_wrapper(x, t, h(h_i), h_i, species_properties, scen_properties);
    end
end
