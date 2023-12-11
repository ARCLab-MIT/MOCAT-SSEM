function [Lambdadot] = launch_func_constant(t, h, species_properties, scen_properties)
    % Adds constant launch rate from species_properties.lambda_constant
    %   t is time from scenario start in years
    %   h is the set of altitudes of the scenario above ellipsoid in km of shell lower edges.
    %   species_properties is a structure with properties for the species
    %   scen_properties.launch_filepath = relative or absolute path to
    %   launch csv file with column for time (year) and one or more launch
    %   rates
    %   scen_properties is a structure with properties for the scenario
    %   Lambdadot is the rate of change in the species in each shell at the 
    %   specified time due to launch If only one value is applied, it is assumed to be true for all shells.
    if length(h) ~= scen_properties.N_shell
        error("Constant launch rate must be specified per altitude shell.")
    end
    
    Lambdadot = ones(scen_properties.N_shell, 1, "sym") .* species_properties.lambda_constant;
    Lambdadot = sym2cell(Lambdadot);

    end