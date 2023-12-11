function eq_index = get_species_index_from_name(scen_properties, species_names)
    % Returns the numerical index for a particular species in
    % scen_properties.species matching to species_name
    % scen_properties is a scen_properties_class with a species array
    % species_name is a string corresponding to the
    % species.species_properties.sym_name property. Does not support
    % indexing for a particular shell.
    eq_index = 0;
    for testi = 1:length(scen_properties.species)
        if scen_properties.species(testi).species_properties.sym_name == species_names
            eq_index = testi;
            break
        end
    end
    if eq_index == 0
        error("Equation index not found for " + species_names)
    end