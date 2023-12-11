classdef indicator_var_class < handle
    properties
        name
        ind_type
        species
        eqs
        indicator_idxs
    end % end properties
    methods
        function obj = indicator_var_class(name, ind_type, species, eqs)
            obj.name = name;
            obj.ind_type = ind_type;
            obj.species = species;
            if exist('eqs', 'var')
                obj.eqs = eqs;
            end
        end
    end %end methods
end %end class
