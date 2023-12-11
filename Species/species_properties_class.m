classdef species_properties_class < handle
    properties
       sym_name
       % Variables
       sym

       % Physical Properties
       Cd
       mass % kg
       mass_lb = 0.00141372 %lower bound of mass class for object binning (inclusive), 1 cm diameter Al sphere
       mass_ub = 100000 %upper bound of mass class for object binning (exclusive, except for top mass bin, which is inclusive)
       radius % m
       A % m^2
       amr % m^2/kg
       beta % m^2/kg
       B 
       density_filepath %For drag

       % Orbit Properties
       slotted = 0 % bool
       slotting_effectiveness = 0 % double [0, 1] where 0 is no col reduction, 1 is perfect col reduction
       disposal_altitude = 0 % km.  If non-zero, satellites that conduct successful PMD will be deposited in circular orbits at this altitude rather than removed from the simulation. Note that this includes satellites at lower altitudes.

       % Capabilities
       drag_effected % bool
       active % bool
       maneuverable % bool
       trackable % bool
       deltat % years
       Pm % double [0, 1], 0 = no Pmd, 1 = full Pm
       alpha % failure rate of collision avoidance vs inactive trackable objects [0 = perfect, 1 = none]
       alpha_active % failure rate of collision avoidance vs active maneuverable objects [0 = perfect, 1 = none]
       RBflag % bool 1 = rocket body, 0 = not rocket body

       % Orbit Raising (Not implemented)
       orbit_raising  % Bool
       insert_altitude % altitude -> semimajor axis
       onstation_altitude
       epsilon % lt_mag
       e_mean

       lambda_multiplier %only for launch_func_fixed_multiplier
       %TODO: Note that lambda_funs, lambda_constant needs to be a cell. Add error checking
       lambda_funs %only for launch_func_gauss
       lambda_constant %only for launch_func_constant
       lambda_python_args

       % For Derelicts
       pmd_linked_species
       pmd_linked_multiplier % Used when multiple shells dispose to one orbit.
    
       % References
       eq_idxs %indexes of species states in X

       % For knowing when to recalculate custom launch functions.
       last_calc_x = nan
       last_calc_t = nan
       last_lambda_dot = nan

       % For looped model
       saved_model_path
       t_plan_max = -1
       t_plan_period = 1
       prev_prop_results

    end % end properties
    methods
        function obj = species_properties_class(struct)
            for field = fieldnames(struct)'    %enumerat fields
                try
                  obj.(field{1}) = struct.(field{1});   %and copy
                catch
                  warning('Error copying field %s', field{1});
                end
            end

        end
    end %end methods
end %end class
