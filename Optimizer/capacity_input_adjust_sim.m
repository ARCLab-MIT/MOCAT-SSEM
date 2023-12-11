function simulation = capacity_input_adjust_sim(simulation, lambda)

    % Takes a simulation and assigns a given column vector for constant lambda 
    % to the S variable.
    % This is used in the wrap_sim part of the optimizer, along with null
    % launch rates to avoid re-building the variables.
    lambda = lambda.'; %Because optimizer feeds in a row, rather than a column
    S_list = get_species_name_list_from_species_cell(simulation, "S");
    S_list = get_species_list_from_names(simulation, S_list);
    S_eq_idxs = []; 
    for i = 1:length(S_list)
        S_eq_idxs = [S_eq_idxs S_list(i).species_properties.eq_idxs];
    end
    
    Su_list = get_species_name_list_from_species_cell(simulation, "Su");
    Su_list = get_species_list_from_names(simulation, Su_list);
    Su_eq_idxs = []; 
    for i = 1:length(Su_list)
        Su_eq_idxs = [Su_eq_idxs Su_list(i).species_properties.eq_idxs];
    end
    
    % Take lambda and add it to the right rows of a copy of the base equations.
    % Then re-write the actual matlabFunction accordingly.
    merged_eqs = zeros(numel(simulation.xdot_eqs), 1, "sym");
    patched_eqs_S = simulation.xdot_eqs(S_eq_idxs) + ones(numel(S_eq_idxs), 1, "sym") .* lambda(1:numel(S_eq_idxs));
    patched_eqs_Su = simulation.xdot_eqs(Su_eq_idxs) + ones(numel(Su_eq_idxs), 1, "sym") .* lambda(numel(S_eq_idxs)+1:end);
    merged_eqs(S_eq_idxs) = patched_eqs_S;
    merged_eqs(Su_eq_idxs) = patched_eqs_Su;
    unchanged_idxs = setdiff(1:length(simulation.xdot_eqs),[S_eq_idxs Su_eq_idxs]);
    merged_eqs(unchanged_idxs) = simulation.xdot_eqs(unchanged_idxs);
    simulation.xdot_fun = matlabFunction(merged_eqs,'Vars',{simulation.scen_properties.var_col});