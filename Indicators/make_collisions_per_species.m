function indicator_var_array = make_collisions_per_species(scen_properties)
    % This function returns an indicator_var_array with a list of indicator
    % vars corresponding to the collisions per year for each species.
    
    dummy_obj = {};
    dummy_obj.scen_properties = scen_properties;
    
    all_col_indicators = indicator_var_class.empty();
    
    for species_i = 1:length(scen_properties.species)
        species_name = scen_properties.species(species_i).species_properties.sym_name;
        spec_col_indicators = indicator_var_class.empty();
        %disp("Loop for " + species_name)
        for pair_i = 1:length(scen_properties.species_pairs) 
            species_1_name = scen_properties.species_pairs(pair_i).species1.species_properties.sym_name;
            species_2_name = scen_properties.species_pairs(pair_i).species2.species_properties.sym_name;
            if species_name == scen_properties.species_pairs(pair_i).species1.species_properties.sym_name || species_name == scen_properties.species_pairs(pair_i).species2.species_properties.sym_name
                ind_name = "collisions_" + species_1_name + "_" + species_2_name;
                spec_pair = [species_1_name species_2_name];
                %disp("     Inner loop for " + spec_pair)
                col_indicator = make_indicator_struct(dummy_obj, ind_name, "collision", spec_pair, []);
                spec_col_indicators = [spec_col_indicators col_indicator];
            end
        end
        ag_col_eqs = sym(zeros(scen_properties.N_shell,1));
        for col_ind_i = 1:length(spec_col_indicators)
            ag_col_eqs = [ag_col_eqs + spec_col_indicators(col_ind_i).eqs];
        end
        spec_ag_col_indc = make_indicator_struct(dummy_obj, species_name + "_aggregate_collisions", "manual", scen_properties.species(species_i), ag_col_eqs);
        all_col_indicators = [all_col_indicators; spec_ag_col_indc];
    end
    
    indicator_var_array = all_col_indicators;