function debris_species = pair_actives_to_debris(scen_properties, active_species, debris_species)
% This function takes a list of species and pairs all active species
% to debris species to properly increment derelict satellites from
% failed PMD and from PMD to a disposal altitude (if this behavior is set
% on a per-species basis for any species using the disposal_altitude param.
% It modifies the species properties for pmd_linked_multiplier and
% pmd_linked_species for debris species.

linked_spec_list = [];
linked_spec_names = [];
for i = 1:numel(active_species)
    cur_species = active_species(i);
    if cur_species.species_properties.active == true
        linked_spec_list = [linked_spec_list cur_species];
        linked_spec_names = [linked_spec_names cur_species.species_properties.sym_name];
    end
end

disp("Pairing the following active species to debris classes for PMD modeling...")
disp(linked_spec_names)

% Assign matching debris increase for a species due to failed PMD
for i = 1:length(linked_spec_list)
    found_mass_match_debris = false; 
    spec_mass = linked_spec_list(i).species_properties.mass;
    for species_i = 1:length(debris_species)
        if spec_mass == debris_species(species_i).species_properties.mass
            debris_species(species_i).pmd_func = @pmd_func_derelict;
            debris_species(species_i).species_properties.pmd_linked_species = [debris_species(species_i).species_properties.pmd_linked_species; linked_spec_list(i)];
            disp("Matched species " ...
                + linked_spec_list(i).species_properties.sym_name ...
                + " to debris species " ...
                +  debris_species(species_i).species_properties.sym_name + ".")
            found_mass_match_debris = true;
        end
    end
    if ~found_mass_match_debris
        warning("No matching mass debris species found for species " + ...
                linked_spec_list(i).species_properties.sym_name + ...
                " with mass " + string(spec_mass) + "." )
    end
end

% Successful PMD to non-zero altitude.
for species_i = 1:length(debris_species)
    cur_deb_spec = debris_species(species_i);

    linked_spec_names = [];
    for i = 1:numel(cur_deb_spec.species_properties.pmd_linked_species)
        cur_species = cur_deb_spec.species_properties.pmd_linked_species(i);
        linked_spec_names = [linked_spec_names cur_species.species_properties.sym_name];
    end
    disp("    Name: " + cur_deb_spec.species_properties.sym_name)
    disp("    pmd_linked_species: " + linked_spec_names)
    %disp("    pmd_linked_multiplier: ")
    %disp("    " + string(cur_deb_spec.species_properties.pmd_linked_multiplier))

    cur_deb_spec_spec_prop = cur_deb_spec. species_properties;
    % num_linked_species =  numel(cur_deb_spec.species_properties.pmd_linked_species);
    % cur_deb_spec_spec_prop.pmd_linked_multiplier = ones(scen_properties.N_shell, num_linked_species);
    % for i = 1:numel(cur_deb_spec_spec_prop.pmd_linked_species)
    %         linked_spec_list = cur_deb_spec_spec_prop.pmd_linked_species;
    %         % Check for circular PMD
    %         % if linked_spec_list(i).species_properties.disposal_altitude ~= 0
    %         %     disp(" Species " + linked_spec_list(i).species_properties.sym_name + ...
    %         %          " has disposal altitude set at " + ...
    %         %          string(linked_spec_list(i).species_properties.disposal_altitude) + ...
    %         %          " km. Adjusting pmd_linked_multiplier for " + ...
    %         %          cur_deb_spec_spec_prop.sym_name +  " accordingly.")
    %         %     disposal_bin = find_alt_bin(linked_spec_list(i).species_properties.disposal_altitude, scen_properties);
    %         %     % Shells above the disposal altitude all dispose to
    %         %     % disposal altitude
    %         %     cur_deb_spec_spec_prop.pmd_linked_multiplier(:,1) = ones(scen_properties.N_shell, 1);
    %         %     cur_deb_spec_spec_prop.pmd_linked_multiplier(disposal_bin, i) = cur_deb_spec_spec_prop.pmd_linked_multiplier(disposal_bin) + sum(cur_deb_spec_spec_prop.pmd_linked_multiplier(disposal_bin + 1:end));
    %         %     cur_deb_spec_spec_prop.pmd_linked_multiplier(disposal_bin + 1:end, i) = 0;
    %         %     disp("Disposal multipliers: " + string(cur_deb_spec_spec_prop.pmd_linked_multiplier))
    %         % end
    % end % Loop over pmd_linked_species
end % Loop over species
