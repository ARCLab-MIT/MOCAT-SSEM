function indicator_var = make_active_loss_per_shell_counter(scen_properties, percentage, per_species)
    % Calculates indicator variable for number of active spacecraft lost in
    % each orbital shell to collision events in a given year.
    % percentage = false uses absolute number, = true gives percentage
    % lost.
    % per_species = true, returns values for each species as independent
    % indicator. per_speices = falses, sums per shell.

    dummy_obj = {};
    dummy_obj.scen_properties = scen_properties;
    species_pairs_classes = scen_properties.species_pairs;
    all_col_indicators = indicator_var_class.empty();

    for species_i = 1:length(scen_properties.species) %Loop for each species pair
        species_1_name = scen_properties.species(species_i).species_properties.sym_name;
        spec_col_indicators = indicator_var_class.empty();
        for pair_i = 1:length(species_pairs_classes) 
            if species_1_name == species_pairs_classes(pair_i).species1.species_properties.sym_name || species_1_name == species_pairs_classes(pair_i).species2.species_properties.sym_name
                species_2_name = species_pairs_classes(pair_i).species2.species_properties.sym_name;
                ind_name = "collisions_" + species_1_name + "_" + species_2_name;
                spec_pair = [species_pairs_classes(pair_i).species1.species_properties.sym_name species_pairs_classes(pair_i).species2.species_properties.sym_name];
                col_indicator = make_indicator_struct(dummy_obj, ind_name, "collision", spec_pair, []);
                spec_col_indicators = [spec_col_indicators col_indicator];
            end
        end
        ag_col_eqs = sym(zeros(scen_properties.N_shell,1));
        for col_ind_i = 1:length(spec_col_indicators)
            ag_col_eqs = [ag_col_eqs + spec_col_indicators(col_ind_i).eqs];
        end
        spec_ag_col_indc = make_indicator_struct(dummy_obj, species_1_name + "_aggregate_collisions", "manual", scen_properties.species(species_i), ag_col_eqs);
        all_col_indicators = [all_col_indicators spec_ag_col_indc];
    end
    
    if per_species == true
        if percentage == false
            indicator_var = all_col_indicators;
        end
        if percentage == true
            for species_i = 1:length(scen_properties.species)
                species_1_name = scen_properties.species(species_i).species_properties.sym_name;
                species_1_totals = scen_properties.species(species_i).species_properties.sym;
                all_col_indicators(species_i).eqs = 100 * all_col_indicators(species_i).eqs./species_1_totals;
            end
        end
    elseif per_species == false
        % Get the active collisions
        ag_active_col_eqs = sym(zeros(scen_properties.N_shell,1));
        for species_i = 1:length(scen_properties.species)
            if scen_properties.species(species_i).species_properties.active
                %disp(scen_properties.species(species_i).species_properties.sym_name)
                ag_active_col_eqs = ag_active_col_eqs + all_col_indicators(species_i).eqs;
            end
        end
    
        if percentage == false
            indicator_var = make_indicator_struct(dummy_obj, "active_aggregate_collisions", "manual", scen_properties.species(species_i), ag_active_col_eqs);
        end
    
        % Get the total active spacecraft
        ag_active_sat_totals = sym(zeros(scen_properties.N_shell,1));
        if percentage == true
            for species_i = 1:length(scen_properties.species)
                if scen_properties.species(species_i).species_properties.active
                    ag_active_sat_totals = ag_active_sat_totals + scen_properties.species(species_i).species_properties.sym;
                end
            end
            perc_eqs = 100*ag_active_col_eqs./ag_active_sat_totals;
            indicator_var = make_indicator_struct(dummy_obj, "active_aggregate_collisions_percentage", "manual", scen_properties.species(species_i), perc_eqs);
        end
    end
