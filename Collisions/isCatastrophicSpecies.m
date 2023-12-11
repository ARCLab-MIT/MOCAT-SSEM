function [overThreshold] = isCatastrophicSpecies(species1, species2, scenario_properties)
    mass1 = species1.species_properties.mass;
    mass2 = species2.species_properties.mass;
    vels = scenario_properties.v_imp2;
    overThreshold = isCatastrophic(mass1, mass2, vels);