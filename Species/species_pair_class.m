classdef species_pair_class
    properties
        name
        species1
        species2
        sigma
        phi
        catastrophic
        nf
        delta
        gammas
        source_sinks
        eqs

    end % end properties
    methods
        function obj = species_pair_class(species1, species2, ...
                                          gammas, source_sinks, ...
                                          scen_properties)
            % This makes the species pair class associated with a collision
            % between species1 and species2. It will then create equations
            % for the collision probability modifiers in gamma and the
            % species in source_sinks.  If the symbolic argument "n_f" is
            % passed, it will be replaced with the n_f value for a
            % collision involving species1 and species2 at each dv in 
            % scen_properties.v_imp2.
            % Inputs:
            %   species1 = species class corresponding to a collision pair
            %   species2 = species class corresponding to a collision pair
            %   gammas = collision rate modifier from things like collision
            %   avoidance, slotting, etc. Generally either a scalar or a N
            %   x M matrix, where N is the number of altitude bins and M is
            %   the number of species with population addition/subtractions
            %   when this collision type occurs.
            %   source_sinks = 1 x M species array for the species with
            %   populations impacted by this collision pair.
            % Outputs:
            %    species_pair_class object with equations that can be fed
            %    to the integrator for this collision pair.

            
            if width(gammas) ~= width(source_sinks) 
                error("Gammas and source_sinks must be the same length")
            end
            obj.name = strcat("species_pair(", ...
                              species1.species_properties.sym_name,",",...
                              species2.species_properties.sym_name, ")");
            obj.species1 = species1;
            obj.species2 = species2;
    
            meter_to_km = 1/1000;

            %Square of impact parameter
            obj.sigma = (species1.species_properties.radius * meter_to_km + ...
                         species2.species_properties.radius * meter_to_km)^2;

            % Scaling based on v_imp, shell volume, and object radii
            obj.phi = pi * scen_properties.v_imp2./(scen_properties.V* meter_to_km^3) * obj.sigma * 86400*365.25; % time duration fix (to 'per year');

            % Check if collision is catastrophic based on 40 j/g criteria
            % and use correct evolve equation to estimate number of
            % fragments above scenario's characteristic length
            obj.catastrophic = isCatastrophicSpecies(species1, species2, ...
                                                     scen_properties);
            % TODO: Note that col code will need to change for vector nf to
            % use nf(k)
            % TODO: Not for D?
            
            % 
            M1 = obj.species1.species_properties.mass;
            M2 = obj.species2.species_properties.mass;
            LC = scen_properties.LC; %TODO: Change when add LNT? Add mass conservation?
            n_f_catastrophic = @(M1,M2) 0.1*LC^(-1.71)*(M1+M2)^(0.75)* ones(size(scen_properties.v_imp2)); % number of fragments generated during a catastrophic collision (NASA standard break-up model). M is the sum of the mass of the objects colliding in kg
            n_f_damaging = @(M1,M2) 0.1*LC^(-1.71)*(min(M1,M2)*scen_properties.v_imp2.^2).^(0.75); % number of fragments generated during a non-catastrophic collision (improved NASA standard break-up model: takes into account the kinetic energy). M is the mass of the less massive object colliding in kg
            if obj.catastrophic
                obj.nf = n_f_catastrophic(M1,M2).';
            else
                obj.nf = n_f_damaging(M1,M2).';
            end
            
            % For each gamma (which modifies collisions for things like
            % collision avoidance, or fragmentation into derelicts and
            % debris, we increment a source or sink species equation.
            obj.gammas = gammas;
            obj.source_sinks = source_sinks;
            eqs = zeros(scen_properties.N_shell, width(scen_properties.species), 'sym'); %Rows are heights, cols are species.
            for i = 1:width(gammas)
                gamma = gammas(:, i);
                % This section finds the index for the species included in
                % source_sink list at a particular index so we know which
                % equation to index. 
                % TODO: Eventually, replace this with a dictionary.

                eq_index = 0;
                for testi = 1:length(scen_properties.species)
                    if scen_properties.species(testi).species_properties.sym_name == source_sinks(i).species_properties.sym_name
                        eq_index = testi;
                        break
                    end
                end
                if eq_index == 0
                    error("Equation index not found for" + source_sinks(i).species_properties.sym_name)
                end

                n_f = sym('n_f', [scen_properties.N_shell, 1]);
                %obj.phi = sym('phi'); %Remove after test
                try
                    eq = gamma .* obj.phi.' .* species1.species_properties.sym .* species2.species_properties.sym;
                catch
                    disp('why')
                end

                eq = subs(eq, n_f, obj.nf);
   
                eqs(:,eq_index) = eqs(:,eq_index) + eq;
            end
            obj.eqs = eqs;
        end
    end %end methods
end %end class