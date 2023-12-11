function species_list = get_species_list_from_names(my_sim, species_names)
    % This utility function takes a simulation_class object and a string
    % list of species names. It returns a list of species classes
    % corresponding to the classes with these names.
    species_list = species.empty();
    for species_i = 1:length(species_names)
        species_name = species_names(species_i);
        for i = 1:length(my_sim.species_list)
            sim_species_name = my_sim.species_list(i).species_properties.sym_name;
            %disp(["sim_species_name:", sim_species_name])
            if strcmp(species_name, sim_species_name)
                species_list(length(species_list)+1) = my_sim.species_list(i);
            end
        end
    end
end