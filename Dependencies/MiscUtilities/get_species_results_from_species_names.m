function spec_list_results = get_species_results_from_species_names(my_sim, species_names)
    % Returns the propagated results for a list of species in
    % scen_properties.species matching to species_name
    % scen_properties is a scen_properties_class with a species array
    % species_name is a string corresponding to the
    % species.species_properties.sym_name property. Does not support
    % indexing for a particular shell.
    idxs = [];
    for i = 1:numel(species_names)
        idx = get_species_index_from_name(my_sim.scen_properties,species_names(i));
        idxs = [idxs idx];
    end

    spec_list_results = my_sim.results.species_array(:,:,idxs);