function spec_name_list = get_species_name_list_from_species_cell(my_sim, species_cell_name)
    % This utility function takes a simulation_class object and a string
    % for a species_cell. It returns a list of species names
    % corresponding to the species in the species_cell.
    spec_name_list = strings(0);
    species_cell = my_sim.scen_properties.species_cell.(species_cell_name);
    for i = 1:length(species_cell)
        spec_name_list(numel(spec_name_list)+1) = species_cell(i).species_properties.sym_name;
    end
end