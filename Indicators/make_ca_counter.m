
function ind_struct = make_ca_counter(scen_properties, primary_species_list, secondary_species_list, varargin)

    % This method makes an indicator variable structure corresponding to
    % the number of collision avoidance maneuvers performed in a given year
    % by the primary species.  It assumes that collision avoidance burden
    % is divided evenly if two active species have a conjunction. Slotting
    % effectiveness is not considered, but should be handled through
    % inclusion in gamma.
    
    % Inputs:
    % scen_properties is a scenario properties object for the model in
    % question
    % primary_species_list is a list of the species for which collision avoidance
    % maneuvers should be counted. Can also be passed as a name of a boolean attribute.
    % If so, will select all species where that attribute is true. "all"
    % will use all species.
    % secondary_species_list is the list of species against which collision
    % avoidance maneuvers should be counted. Can also be passed as a name of a boolean attribute.
    % If so, will select all species where that attribute is true. "all"
    % will use all species.
    % per_species is a bool (default false) for whether maneuver counts should be returned
    % in aggregate for all primary species, or on a per-species basis.
    % per_spacecraft is a bool (default false) for whether maneuver counts
    % should be returned per spacecraft instead of per species in each
    % altitude bin.
    % ind_name is a string that overrides default naming options.



    % Warning, assumes gamma for both pairs in collision will be symmetric.
    % This is true in our sample models, but could be different if one
    % assumed that the probability of loss was different across species.

    default_per_species = false;
    default_per_spacecraft = false;
    default_ind_name = "";
    p = inputParser;
    addRequired(p, 'scen_properties', @(x) class(x) == "scen_properties_class")
    addRequired(p, 'primary_species_list')
    addRequired(p, 'secondary_species_list')
    addOptional(p,'per_species',default_per_species, @islogical);
    addOptional(p,'ind_name', default_ind_name);
    addOptional(p, 'per_spacecraft', default_per_spacecraft, @islogical)
    parse(p, scen_properties, primary_species_list, secondary_species_list, varargin{:})
    per_species = p.Results.per_species;
    per_spacecraft = p.Results.per_spacecraft;
    ind_name = p.Results.ind_name;
    primary_species_list_attribute = "";
    secondary_species_list_attribute = "";

    if isstring(primary_species_list)
        primary_species_list_attribute = primary_species_list;
        primary_species_list = [];
        for species_i = 1:length(scen_properties.species)
            if primary_species_list_attribute == "all"
                primary_species_list = scen_properties.species;
                break
            end
            if scen_properties.species(species_i).species_properties.(primary_species_list_attribute)
                %disp(scen_properties.species(species_i).species_properties.sym_name)
                primary_species_list = [primary_species_list scen_properties.species(species_i)];
            end
        end
    end

    if isstring(secondary_species_list)
        secondary_species_list_attribute = secondary_species_list;
        secondary_species_list = [];
        for species_i = 1:length(scen_properties.species)
            if secondary_species_list_attribute == "all"
                secondary_species_list = scen_properties.species;
                break
            end
            if scen_properties.species(species_i).species_properties.(secondary_species_list_attribute)
                %disp(scen_properties.species(species_i).species_properties.sym_name)
                secondary_species_list = [secondary_species_list scen_properties.species(species_i)];
            end
        end
    end

    % Get species name cell arrays from the species lists
    primary_species_names = [];
    secondary_species_names = [];
    for name_i = 1:length(primary_species_list)
        primary_species_names = [primary_species_names primary_species_list(name_i).species_properties.sym_name];
    end
    for name_i = 1:length(secondary_species_list)
        secondary_species_names = [secondary_species_names secondary_species_list(name_i).species_properties.sym_name];
    end
    if  ind_name == ""
        ind_name = "col_avoidance_maneuvers_pri_";
        if primary_species_list_attribute ~= ""
            ind_name = ind_name + primary_species_list_attribute ;
        else
            ind_name = ind_name + strjoin(primary_species_names,'_') ;
        end
        if secondary_species_list_attribute ~= ""
            ind_name = ind_name + "_sec_" + secondary_species_list_attribute ;
        else
            ind_name = ind_name + "_sec_" + strjoin(secondary_species_names,'_') ;
        end
        if per_spacecraft
            ind_name = ind_name + "_per_spacecraft";
        end
    end
    
    ind_eqs = {};
    full_spec_list = [primary_species_list secondary_species_list];
    for species_i = 1:length([primary_species_list secondary_species_list])
        species = full_spec_list(species_i);
        species_name = species.species_properties.sym_name;
        if species.species_properties.maneuverable && species.species_properties.active
            ind_eqs.(species_name) = sym(zeros(scen_properties.N_shell,1));
        end
    end

    % Filter the col pair list
    for species_i = 1:length(primary_species_list)
        primary_species = primary_species_list(species_i);
        primary_species_name = primary_species.species_properties.sym_name;
        %disp("Checking species: " + primary_species_name)
        for species_pair_class_index = 1:length(scen_properties.species_pairs)
            pair_primary_name = scen_properties.species_pairs(species_pair_class_index).species1.species_properties.sym_name;
            pair_secondary_name = scen_properties.species_pairs(species_pair_class_index).species2.species_properties.sym_name;
            %disp("    Checking pair: " + pair_primary_name + " " + pair_secondary_name)
            s1_prim = primary_species_name == pair_primary_name;
            s2_prim = primary_species_name == pair_secondary_name;
            if s1_prim
                sec_in_sec_list = any(strcmp(secondary_species_names,pair_secondary_name));
                secondary_species_name = pair_secondary_name;
                %disp("        S1 is primary")
            end
            if s2_prim
                sec_in_sec_list = any(strcmp(secondary_species_names,pair_primary_name));
                secondary_species_name = pair_primary_name;
                %disp("        S2 is primary")
            end
            primary_and_sec_in_pair = (s1_prim || s2_prim) && (sec_in_sec_list);
            if primary_and_sec_in_pair
                %disp("        S1 and S2 in primary or secondary list")
                pair = scen_properties.species_pairs(species_pair_class_index);
                %disp(pair_primary_name + " " + pair_secondary_name)
                %disp(pair)
                % Non-Gamma
                intrinsic_collisions = pair.phi.' .* pair.species1.species_properties.sym .* pair.species2.species_properties.sym;
                
                both_man = pair.species1.species_properties.maneuverable && pair.species2.species_properties.maneuverable;
                both_act = pair.species1.species_properties.active && pair.species2.species_properties.active;
                s1_man_act = pair.species1.species_properties.maneuverable && pair.species1.species_properties.active;
                s2_man_act = pair.species2.species_properties.maneuverable && pair.species2.species_properties.active;
                s1_trackable = pair.species1.species_properties.trackable;
                s2_trackable = pair.species2.species_properties.trackable;

                one_man_act = xor(s1_man_act, s2_man_act);
                
                if  both_man && both_act % Avoid double-counting maneuvers for active on active.
                    %disp("            Both species maneuver and active.")
                    if pair.species1.species_properties.slotted && pair.species2.species_properties.slotted
                        %TODO: Check This.
                        %disp("                Both species slotted.")
                        slotting_effectiveness = min(pair.species1.species_properties.slotting_effectiveness, pair.species2.species_properties.slotting_effectiveness);
                        ind_eqs.(primary_species_name)    = ind_eqs.(primary_species_name) + 0.5* (1+pair.gammas(:,1)) * (1/slotting_effectiveness) .* intrinsic_collisions; 
                        ind_eqs.(secondary_species_name) =  ind_eqs.(secondary_species_name) + 0.5* (1+pair.gammas(:,1)) * (1/slotting_effectiveness).* intrinsic_collisions;
                    else
                        ind_eqs.(primary_species_name)    = ind_eqs.(primary_species_name) + 0.5* (1+pair.gammas(:,1)) .* intrinsic_collisions; 
                        ind_eqs.(secondary_species_name) =  ind_eqs.(secondary_species_name) + 0.5* (1+pair.gammas(:,1)) .* intrinsic_collisions;
                    end

                elseif one_man_act
                    % Make Eqs
                    if s1_man_act && s2_trackable
                        %disp("            Only s1 is maneuverable and active. S2 is trackable.")
                        ind_eqs.(pair_primary_name) = ind_eqs.(pair_primary_name) + (1+pair.gammas(:,1)) .* intrinsic_collisions;
                    elseif s2_man_act && s1_trackable
                        %disp("            Only s2 is maneuverable and active. S1 is trackable.")
                        ind_eqs.(pair_secondary_name) = ind_eqs.(pair_secondary_name) + (1+pair.gammas(:,1)) .* intrinsic_collisions;
                    end
                % Make gamma eqs.
                end
            end
        end
    end

    if per_species == false
        ag_man_counts = sym(zeros(scen_properties.N_shell,1));
        fn = fieldnames(ind_eqs);
        for k=1:numel(fn)
            ag_man_counts = ag_man_counts + ind_eqs.(fn{k});
        end
        if per_spacecraft == true
            ag_man_sat_totals = sym(zeros(scen_properties.N_shell,1));
            for species_i = 1:length(scen_properties.species)
                if scen_properties.species(species_i).species_properties.maneuverable
                    ag_man_sat_totals = ag_man_sat_totals + scen_properties.species(species_i).species_properties.sym;
                end
            end
            ag_man_counts = ag_man_counts./ag_man_sat_totals;
        end
        ind_struct = indicator_var_class(ind_name, "manual", [], ag_man_counts);
    elseif per_species == true
        ind_struct = indicator_var_class.empty();
        fn = fieldnames(ind_eqs);
        for k=1:numel(fn)
            if per_spacecraft == true
                spec_index = get_species_index_from_name(scen_properties, fn{k});
                spec_man_indc = indicator_var_class(fn{k} + "_maneuvers_per_spacecraft", "manual", [], ind_eqs.(fn{k})./scen_properties.species(spec_index).species_properties.sym);
            elseif per_spacecraft == false
                spec_man_indc = indicator_var_class(fn{k} + "_maneuvers", "manual", [], ind_eqs.(fn{k}));
            end
            ind_struct = [ind_struct; spec_man_indc];
        end
    end

