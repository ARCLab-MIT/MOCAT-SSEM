function [Cpmddot] = pmd_func_none(t, h, species_properties, scen_properties)
    %   PMD function for placeholder classes with no PMD, e.g. debris
    %   t is time from scenario start in years (unused)
    %   h is the height above ellipsoid in km (unused)
    %   species_properties is a structure with properties for this species
    %   Cpmdot is the rate of change in the species due to post-mission
    %   disposal, an N_shell x 1 matrix.
    
    Cpmddot = zeros(scen_properties.N_shell, 1, "sym") .* species_properties.sym(:); %Species is not strictly necessary

    end