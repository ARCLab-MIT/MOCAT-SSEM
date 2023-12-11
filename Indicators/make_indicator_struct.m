function new_indicator_var_class = make_indicator_struct(obj, name, ind_type, species, eqs)
    % Helper method that creates a indicator variable structure. This is
    % intended to be called during the simulation_class build_model method
    % for each indicator variable intended to be added to the simulation.
    %   obj is a my_sim object.
    %
    %   name is the name for the indicator variable
    %   species is a string array corresponding to the species for the 
    %   indicator variables.
    % 
    %   TODO: smartly parse species from list rather than require exact
    %   match
    %   species must be passed, but is ignored if type is manualmy_sim.species_pairs
    %   Three options are supported
    %   First, the code checks to see if the name of the variable matches
    %   to the sym_name of a species, e.g. "S". If a species is split into multiple mass sub-species,
    %   This also will include and add values for all subspecies.
    %   Second. If there is not a match to a general species, the code
    %   checks if there is an exact match to a subspecies.
    %   Third, the code checks if there is a match to the name of a
    %   particular variable associated with a single shell of a particular
    %   species, e.g. "Su_6".
    %
    %   eqs is a symbolic value corresponding to the equations to use for a
    %   particular constraint. This is an optional argument, intended to be
    %   passed only for ind_type = "manual"

    %Test cases: 
    % 1) pass type = manual and no eqs
    % 2) pass type = "drag" and eqs

    % Note that collision and mitigated conjunctions use first gamma from
    % species pair
    valid_types = ["collision", "successful PMD", "failed PMD", "mitigated conjunctions", "manual"];

    % Error checking for manual or not
    if strcmp(ind_type, "manual") && ~exist('eqs','var')
        error("Argument eqs must be passed for indicator variable with type ''manual''.")
    end

    if ~strcmp(ind_type, "manual") && exist('manual_eqs','var')
        warning("Argument eqs passed for indicator with non-manual type. Overriding equations to use passed value. To remove this error, use type = ''manual''.")
    end

    % Check if ind_type is a supported type

    if ~ismember(ind_type,valid_types)
        error("Passed type " + string(ind_type) + " is not a supported type. Please use values from" + valid_types + ".")
    end

    % Build the equations
    ind_struct.name = name;
    ind_struct.type = ind_type;
    ind_struct.species = species;
    
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
                break
            end   
        end
        if ~exist('pair','var')
            error("No matching species pair has been found for: ''" + strjoin(ind_struct.species, "'', ''") + "'' .")
        end

        % Non-Gamma
        intrinsic_collisions = pair.phi.' .* pair.species1.species_properties.sym .* pair.species2.species_properties.sym;
        
        % Make Eqs
        ind_struct.eqs = -1 * pair.gammas(:,1) .* intrinsic_collisions; % Negative 1 is to counteract decrease to quanity to provide positive number of collisions

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
            error("No matching species pair has been found for: ''" + strjoin(ind_struct.species, "'', ''") + "'' .")
        end

        % Non-Gamma
        intrinsic_collisions = pair.phi.' .* pair.species1.species_properties.sym .* pair.species2.species_properties.sym;
        
        % Make Eqs
        ind_struct = (1-pair.gammas(:,1)) .* intrinsic_collisions;
%     elseif strcomp(ind_type, "successful PMD")
% 
% 
%         obj.species_list(i).pmd_func(t, obj.scen_properties.HMid, obj.species_list(i).species_properties, obj.scen_properties)
%     elseif strcomp(ind_type, "failed PMD")

    elseif strcmp(ind_type, "manual")
        ind_struct.eqs = eqs;
    end
    new_indicator_var_class = indicator_var_class(ind_struct.name, ind_struct.type, ind_struct.species, ind_struct.eqs);
end