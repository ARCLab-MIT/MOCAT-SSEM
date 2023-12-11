function ind_struct = make_indicator_eqs(obj, ind_struct)
    % Helper method that creates the indicator_var_class.eqs for non-manual
    % types and and does some validation. This is intended to be called 
    % during the simulation_class build_model method for each 
    % indicator_var_class intended to be added to the simulation.

    % TODO: Add error catching to indicate no support for time-varying
    % equations.indicator_var_list
    name = ind_struct.name;
    ind_type = ind_struct.ind_type;
    species = ind_struct.species;
    if isprop(ind_struct, 'eqs')
        eqs = ind_struct.eqs;
    end

    % Note that collision and mitigated conjunctions use first gamma from
    % species pair

    valid_types = ["collision", "successful PMD", "failed PMD", "mitigated conjunctions", "manual"];

    % Error checking for manual or not
    if strcmp(ind_type, "manual") && (~isprop(ind_struct, eqs) || ~isempty(ind_struct.eqs))
        error("Argument eqs must be passed for indicator variable with type ''manual''.")
    end

    if ~strcmp(ind_type, "manual") && isprop(ind_struct, eqs) && ~isempty(ind_struct.eqs)
        warning("Argument eqs passed for indicator with non-manual type. Overriding equations to use passed value. To remove this warning, use type = ''manual''.")
    end

    % Check if ind_type is a supported type

    if ~ismember(ind_type,valid_types)
        error("Passed type " + string(ind_type) + " is not a supported type. Please use values from" + valid_types + ".")
    end

    % Build the equations

    if strcmp(ind_type, "collision")
        if length(species) ~= 2
            error("Exactly two species must be provided for collision indicator variables.")
        end
        
        % Find the correct pair
        for pair_i = 1:length(obj.scen_properties.species_pairs)
            test_pair = obj.scen_properties.species_pairs(pair_i);
            species1match = strcmp(test_pair.species1.species_properties.sym_name, ind_struct.species(1));
            species2match = strcmp(test_pair.species2.species_properties.sym_name, ind_struct.species(2));
            if species1match && species2match
                pair = test_pair;
            end   
        end
        if ~exist('pair','var')
            error("No matching species pair has been found for: ''" + strjoin(ind_struct.species, "'', ''") + "'' .")
        end

        % Non-Gamma
        intrinsic_collisions = pair.phi.' .* pair.species1.species_properties.sym .* pair.species2.species_properties.sym;
        
        % Make Eqs
        ind_struct.eqs = pair.gammas(:,1) .* intrinsic_collisions;

    elseif strcmp(ind_type, "mitigated conjunction")
        if length(species) ~= 2
            error("Exactly two species must be provided for mitigated conjunction indicator variables.")
        end
        
        % Find the correct pair
        for pair_i = 1:length(obj.scen_properties.species_pairs)
            test_pair = obj.scen_properties.species_pairs(pair_i);
            species1match = strcmp(test_pair.species1.species_properties.sym_name, ind_struct.species(1));
            species2match = strcmp(test_pair.species2.species_properties.sym_name, ind_struct.species(2));
            if species1match && species2match
                pair = test_pair;
            end   
        end
        if ~exist('pair','var')
            error("No matching species pair has been found for: ''" + strjoin(ind_var_class.species, "'', ''") + "'' .")
        end

        % Non-Gamma
        intrinsic_collisions = pair.phi.' .* pair.species1.species_properties.sym .* pair.species2.species_properties.sym;
        
        % Make Eqs
        ind_struct.eqs = (1-pair.gammas(:,1)) .* intrinsic_collisions;
    elseif strcmp(ind_type, "successful PMD")
        if length(ind_struct.species)>1
            error("Only single species arguments are currently supported for the ''successful PMD'' type.")
        end
        
        species = get_species_list_from_names(obj, ind_struct.species);
        
        if ~any(ismember(functions(species.pmd_func).function, ['pmd_func', 'pmd_func_sat', 'pmd_func_derelict', 'pmd_func_none']))
            warning("Custom PMD function detected. Warning: only constant (non-time-varying) PMD functions are currently supported for indicator variables.");
        end
        eqs = species.pmd_func(0, obj.scen_properties.HMid, species.species_properties, obj.scen_properties);
        ind_struct.eqs = eqs;
        % obj.species_list(i).pmd_func(t, obj.scen_properties.HMid, obj.species_list(i).species_properties, obj.scen_properties)
    elseif strcmp(ind_type, "failed PMD")
        if length(ind_struct.species)>1
            error("Only single species arguments are currently supported for the ''unsuccessful PMD'' type.")
        end
        species = get_species_list_from_names(obj, ind_struct.species);
        if ~any(ismember(functions(species.pmd_func).function, ['pmd_func', 'pmd_func_sat', 'pmd_func_derelict', 'pmd_func_none']))
            warning("Custom PMD function detected. Warning: only constant (non-time-varying) PMD functions are currently supported for indicator variables.");
        end
        eqs = zeros(obj.scen_properties.N_shell, 1, 'sym');
        for species_i = 1:length(obj.scen_properties.species)
            species_test = obj.scen_properties.species(species_i);
            %disp(["Checking", species_test.species_properties.sym_name, "Linked Species"])
            %disp(species_test.species_properties.pmd_linked_species)
            for species_pmd_i = 1:length(species_test.species_properties.pmd_linked_species) % Per linked species.
                linked_species = species_test.species_properties.pmd_linked_species(species_pmd_i);
                if strcmp(linked_species.species_properties.sym_name, ind_struct.species)
                    % disp("Match found for " + obj.scen_properties.species(species_i).species_properties.sym_name + ".")
                    % If there are multiple pmd_linked_species, we only
                    % want to use the correct one.
                    actualPMD = obj.scen_properties.species(species_i).species_properties.pmd_linked_species;
                    obj.scen_properties.species(species_i).species_properties.pmd_linked_species = linked_species;
                    eqs = eqs + species_test.pmd_func(0, ...
                                                      obj.scen_properties.HMid, ...
                                                      obj.scen_properties.species(species_i).species_properties, ...
                                                      obj.scen_properties);
                    obj.scen_properties.species(species_i).species_properties.pmd_linked_species = actualPMD;
                end
            end
        end  
    elseif strcmp(ind_type, "manual")
        ind_struct.eqs = eqs;
    end

end