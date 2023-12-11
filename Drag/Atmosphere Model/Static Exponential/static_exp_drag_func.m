function [p] = static_exp_drag_func(t, h, species_properties, scen_properties)
% static_exp_drag_func Wrapper for densityexp to be used by species
%   constructor.
%   t is time from scenario tatart in years (unused)
%   h is the height above ellipsoid in km.
%   species_properties is a structure with properties for the species
%   (unused)
%   scen_properties is a structure with properties for the scenario
%   (unused)
%   density p is returned in kg/km^3
p = densityexp(h);
end