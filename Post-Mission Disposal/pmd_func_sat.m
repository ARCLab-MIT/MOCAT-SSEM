function [Cpmddot] = pmd_func_sat(t, h, species_properties, scen_properties)
    %   PMD function for active classes that decay to a derelict. Only
    %   decreases the class based on assumed successful PMD
    %   t is time from scenario start in years (unused)
    %   h is the height above ellipsoid in km (unused)
    %   species_properties is a structure with properties for this species
    %   scen_properties is a structure with properties for the scenario
    %   Cpmdot is the rate of change in the species due to post-mission
    %   disposal, an N_shell x 1 matrix.
    Cpmddot = zeros(scen_properties.N_shell, 1, 'sym');
    for k=1:scen_properties.N_shell
        Cpmddot(k, 1) = (-1/species_properties.deltat) * species_properties.sym(k);
    end
    