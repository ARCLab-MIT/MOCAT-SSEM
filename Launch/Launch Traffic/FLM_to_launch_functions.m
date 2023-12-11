
function scen_properties = FLM_to_launch_functions(FLM_steps, scen_properties)
    % This function uses the FLM created by ADEPT_Traffic_Model and makes a
    % new launch function for each species based on the returned future
    % launch model.  scen_properties is returned for convenience, but the
    % passed object is modified in place as it is a handle class.

    % Take the FLM_steps and make into launch rates for the species.
    for spec_index = 1:length(scen_properties.species)
       spec_names(spec_index) = scen_properties.species(spec_index).species_properties.sym_name;
    end
    
    if length(unique(round(diff(scen_properties.scen_times),5))) == 1
        time_step = unique(round(diff(scen_properties.scen_times),5));
    else
        error("FLM to Launch Function is not stet up for variable time step runs.")
    end

    for ii = 1:length(scen_properties.species)
        sel_species = scen_properties.species(ii);
        species_names = sel_species.species_properties.sym_name;
        %spec_index = get_species_index_from_name(scen_properties, species_names);
    
        spec_FLM = array2table(zeros(scen_properties.N_shell,length(scen_properties.scen_times)), ...
                           RowNames = string(1:scen_properties.N_shell), ...
                           VariableNames = string(scen_properties.scen_times));
        
        % Make the FLM array (Nshell x NTimes)
        % Assume no launches at Tstep 1, as the pop is considered x0
        for jj = 2:length(FLM_steps) % JJ iterates over times
            spec_FLM(:,jj) = FLM_steps{jj}(:,ii); % ii is species
        end
        
        % Scale sats to sats/year and make an interp
        Lambdadot = cell(size(scen_properties.N_shell,2), 0);
        for hh = 1:scen_properties.N_shell
            x = scen_properties.scen_times;
            v = spec_FLM{hh,:}.'/time_step; % Convert from sats to sats/year
            Lambdadot{hh, 1} = @(X, t) interp1(x, v, t, 'linear');
        end
        sel_species.species_properties.lambda_funs = Lambdadot;
        sel_species.launch_func = @launch_func_lambda_fun;
end