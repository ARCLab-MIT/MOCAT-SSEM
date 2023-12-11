function [Lambdadot] = launch_func_lambda_fun(t, h, species_properties, scen_properties)
    % Rename of launch_func_gauss with same function.  Takes lambda from functions stores under species.lambda_funcs as a N_shell x 1 array.
    %   t is time from scenario start in years
    %   h is the altitude above ellipsoid in km of shell lower edges. either scalar or array
    %   species_properties is a structure with properties for the species
    %   scen_properties.launch_filepath = relative or absolute path to
    %   launch csv file with column for time (year) and one or more launch
    %   rates
    %   scen_properties is a structure with properties for the scenario
    %   Lambdadot is the rate of change in the species in each shell at the
    %   specified time due to launch, as cell array.
  
    h_ind = find(scen_properties.HMid==h);
    disp(species_properties.sym_name)
    Lambdadot = species_properties.lambda_funs(h_ind, 1);
   
    end