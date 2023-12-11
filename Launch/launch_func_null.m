function [Lambdadot] = launch_func_null(t, h, species_properties, scen_properties)
    % Takes discrete launch function from file.
    %   t is time from scenario start in years
    %   h is the altitude above ellipsoid in km of shell lower edges. either scalar or array
    %   species_properties is a structure with properties for the species
    %   scen_properties.launch_filepath = relative or absolute path to
    %   launch csv file with column for time (year) and one or more launch
    %   rates
    %   scen_properties is a structure with properties for the scenario
    %   Lambdadot is the rate of change in the species in each shell at the specified time due to launch.
  
    Lambdadot = sym2cell(zeros(scen_properties.N_shell, 1, "sym") .* species_properties.sym(:));
    end