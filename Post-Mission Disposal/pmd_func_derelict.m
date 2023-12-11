function [Cpmddot] = pmd_func_derelict(t, h, species_properties, scen_properties)
    %   PMD function for derelict class. Only
    %   increases the class based on assumed failed PMD from species in
    %   species_properties.linked_species.
    %   t is time from scenario start in years (unused)
    %   h is the height above ellipsoid in km (unused)
    %   species_properties is a structure with properties for this species
    %   scen_properties is a structure with properties for the scenario
    %   Cpmdot is the rate of change in the species due to post-mission
    %   disposal, an N_shell x 1 matrix.
    num_linked_species = length(species_properties.pmd_linked_species);

    Cpmddot = zeros(scen_properties.N_shell, num_linked_species, 'sym');
    for i = 1:num_linked_species
        
        species = species_properties.pmd_linked_species(i);
        Pm = species.species_properties.Pm;
        deltat = species.species_properties.deltat; 
        
        % Matrix, rows are shell vars, and cols are linked species
        Cpmddot(:, i) = (1-Pm)/deltat * species.species_properties.sym(:); % Failed PMD

        % If there's disposal at an altitude instead of removal from
        % scenario
        if species.species_properties.disposal_altitude ~= 0
            disp(" Species " + species.species_properties.sym_name + ...
                     " has disposal altitude set at " + ...
                     string(species.species_properties.disposal_altitude) + ...
                     " km. Adjusting pmd_linked_multiplier for " + ...
                     species_properties.sym_name +  " accordingly.")
            disposal_bin = find_alt_bin(species.species_properties.disposal_altitude, scen_properties);
            % Shells above the disposal altitude all dispose to disposal altitude
            species_properties.pmd_linked_multiplier(:,1) = ones(scen_properties.N_shell, 1);
            species_properties.pmd_linked_multiplier(disposal_bin + 1:end, i) = 0;
            disp("Disposal multipliers: " + string( species_properties.pmd_linked_multiplier))
            species_properties.pmd_linked_multiplier(disposal_bin) = species_properties.pmd_linked_multiplier(disposal_bin);
                            %Failed PMD     and successful PMD below threshold stay
                            %in their shell
            Cpmddot(:, i) = Cpmddot(:, i) + Pm/deltat * species_properties.pmd_linked_multiplier(:,i) .* species.species_properties.sym(:);

            % Above moves to disposal altitude
            Cpmddot(disposal_bin, i) =  Cpmddot(disposal_bin, i) + Pm/deltat * sum(species.species_properties.sym(disposal_bin+1:end));
        end
    end

    Cpmddot = sum(Cpmddot, 2);
  
end