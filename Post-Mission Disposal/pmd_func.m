function [Cpmddot] = pmd_func_satellite(t, h, species_properties, scen_properties)
    % static_exp_drag_func Wrapper for densityexp to be used by species
    %   constructor.
    %   t is time from scenario start in years
    %   h is the height above ellipsoid in km.
    %   species_properties is a structure with properties for the species
    %   (unused)
    %   scen_properties is a structure with properties for the scenario
    %   Cpmdot is the rate of change in the species due to post-mission disposal.
    Cpmddot = (-1/species_properties.deltat) * species_properties.sym(k);
    