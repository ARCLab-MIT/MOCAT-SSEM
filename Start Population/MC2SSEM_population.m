function [popSSEM] = MC2SSEM_population(popMC,VAR)

% The function creats the initial population used by SSEM (divided by
% species and altitude bin) starting from the MC population. The structure
% param contains the parameters that need to be set.
%---
% Authors: Davide Gusmini, Andrea D'Ambrosio, MIT 10/13/2022
% Modified by Miles Lifson to work for MOCAT-SSEM
%---
    
    % species: [S,D,N,Su,B,U];
    ind = find(VAR.species_types~=0);
    % ind

    % Read species, find mass bin edges.
    
    % Read popMC
    n_sats = length(popMC);
    ind_s = 1;   popMC_S = {};
    ind_d = 1;   popMC_D = {};
    ind_n = 1;   popMC_N = {};
    ind_su = 1;  popMC_Su = {};
    ind_b = 1;   popMC_B = {};
    ind_u = 1;   popMC_U = {};
    
    %TODO: re-write so it pulls species list from VAR.species_cell, handles
    %priority order for object_class values depending on which species exist
    
    potential_payload_classes = {'Payload', 'Payload Mission Related Object', 'Other Mission Related Object', 'Rocket Mission Related Object'};
    debris_classes = {'Other Debris', 'Payload Debris',  'Payload Fragmentation Debris', 'Rocket Debris' , 'Rocket Fragmentation Debris'};
    untracked_debris_classes = {'Untracked Debris'};
    rocket_body_classes = {'Rocket Body'};
    
    for k = 1:n_sats
        % Outside boundaries TODO: understand why this is not an || ?
        if (popMC{k}.a * VAR.re - VAR.re) < VAR.h_min && (popMC{k}.a * VAR.re - VAR.re) > VAR.h_max
            continue
        else %classify it
            % Payloads and derelicts
            if any(strcmpi(popMC{k}.objectclass, potential_payload_classes))
                if ~isempty(ind(ind == 1)) && ~isempty(ind(ind == 4)) % S,Su exist
                    % disp("4S")
                    if popMC{k}.controlled == 1 && popMC{k}.constel == 0 % active payload (All assumed unslotted for now?)
                        if ~isempty(ind(ind == 4))
                            popMC_Su{ind_su} = popMC{k};
                            ind_su = ind_su+1;
                        end
                    elseif popMC{k}.controlled == 1 && popMC{k}.constel == 1 % Constellation
                        if ~isempty(ind(ind == 1))
                            popMC_S{ind_s} = popMC{k};
                            ind_s = ind_s+1;
                        end
                    elseif popMC{k}.controlled == 0 % inactive payload (derelict)
                        if ~isempty(ind(ind == 2))
                            popMC_D{ind_d} = popMC{k};
                            ind_d = ind_d+1;
                        end
                    end
                else 
                    % disp("3")
                    if popMC{k}.controlled == 1
                        if ~isempty(ind(ind == 4))
                            popMC_Su{ind_su} = popMC{k};
                            ind_su = ind_su+1;
                        end
                    else
                        if ~isempty(ind(ind == 2))
                            popMC_D{ind_d} = popMC{k};
                            ind_d = ind_d+1;
                        end
                    end
                end
            elseif any(strcmpi(popMC{k}.objectclass, debris_classes))
                if ~isempty(ind(ind == 3))
                    popMC_N{ind_n} = popMC{k};
                    ind_n = ind_n+1;
                end
            elseif any(strcmpi(popMC{k}.objectclass,rocket_body_classes)) % Rocket Bodies
                if ~isempty(ind(ind == 5))
                    popMC_B{ind_b} = popMC{k};
                    ind_b = ind_b+1;
                end
            elseif any(strcmpi(popMC{k}.objectclass,untracked_debris_classes)) % Untracked Debris
                if ~isempty(ind(ind == 6))
                    popMC_U{ind_U} = popMC{k};
                    ind_u = ind_u+1;
                end
            end
        end
    end
        
    % TODO: Redo the species split so that I'm not making U as a subset of
    % N.  It leads to a lot of complicated processing here for no reason. 

    % Binning
    count_list = {};
    pop_list = {};
    speciesName_list = {};
    if isfield(VAR.species_cell, 'S')
        s_count = zeros(VAR.N_shell, length(VAR.species_cell.S)); 
        count_list{length(count_list)+1} = s_count;
        pop_list{length(pop_list)+1} = popMC_S;
        speciesName_list{length(speciesName_list)+1} = "S";
    end
    if isfield(VAR.species_cell, 'D')
        d_count = zeros(VAR.N_shell, length(VAR.species_cell.D));
        count_list{length(count_list)+1} = d_count;
        pop_list{length(pop_list)+1}  = popMC_D;
        speciesName_list{length(speciesName_list)+1}= "D";
    end
    if isfield(VAR.species_cell, 'N')
        n_count = zeros(VAR.N_shell, length(VAR.species_cell.N)); 
        count_list{length(count_list)+1} = n_count;
        pop_list{length(pop_list)+1}  = popMC_N;
        speciesName_list{length(speciesName_list)+1} = "N";
    end
    if isfield(VAR.species_cell, 'Su')
        su_count = zeros(VAR.N_shell, length(VAR.species_cell.Su)); 
        count_list{length(count_list)+1} = su_count;
        pop_list{length(pop_list)+1}  = popMC_Su;
        speciesName_list{length(speciesName_list)+1} = "Su";
    end
    if isfield(VAR.species_cell, 'B')
        b_count = zeros(VAR.N_shell, length(VAR.species_cell.B));
        count_list{length(count_list)+1} = b_count;
        pop_list{length(pop_list)+1}  = popMC_B;
        speciesName_list{length(speciesName_list)+1} = "B";
    end
    if isfield(VAR.species_cell, 'U')
        u_count = zeros(VAR.N_shell, length(VAR.species_cell.U)); 
        count_list{length(count_list)+1} = u_count;
        pop_list{length(pop_list)+1}  = popMC_U;
        speciesName_list{length(speciesName_list)+1} = "U";
    end
    
    % popMC_S 
    % popMC_D 
    % popMC_N 
    % popMC_Su
    % popMC_B 
    % popMC_U 
    % size(count_list)
    
    %% Put the objects into the shell, species, and mass bins based on their info.
    
    for k = 1:VAR.N_shell-1
        for i = 1:length(count_list)
            %ind = count_list{i};
            pop = pop_list{i};
            speciesName = speciesName_list{i};
            for ind = 1:length(pop)
                a_t_shell = pop{ind}.a * VAR.re - VAR.re;
                if a_t_shell>=VAR.R02(k) && a_t_shell<VAR.R02(k+1)
                    species_cell = getfield(VAR.species_cell,speciesName);
                    for mass_i = 1:length(species_cell)
                        if (pop{ind}.mass > species_cell(1,mass_i).species_properties.mass_lb && ...
                            pop{ind}.mass < species_cell(1,mass_i).species_properties.mass_ub )
    %                         disp(["Obj Mass " + num2str(pop{ind}.mass), ...
    %                               " Species mass " + num2str(species_cell(mass_i).species_properties.mass), ...
    %                               " mass lb " + num2str(species_cell(mass_i).species_properties.mass_lb), ...
    %                               " mass ub " + num2str(species_cell(mass_i).species_properties.mass_ub)])
    %                         disp("    Species Found")
                            break %Found right mass bin
                        else
                           continue % Keep looking
                        end
                    end
                    count_list{i}(k, mass_i) = count_list{i}(k, mass_i)+1; % count the number of active satellites in each shell
                end % assign to right shell and species
            end % sat loop
        end % species loop
    end % alt loop
    
    % Assemble list of species in VAR.species (order used by model)
    for i = 1:length(VAR.species)
        name = VAR.species(i).species_properties.sym_name;
        name = strsplit(name, "_");  %TODO this may break if single species are created with _ in the filename
        name = name(1);
        name_order(i) = name;
    end
    
    %Index order
    speciesName_list_idx_list = [];
    for i = 1:length(name_order)
        speciesName_idx = find(strcmp([speciesName_list{:}], name_order(i))); % single line engine;
        % If it's not already in the list of indexes to assemble count_list, add
        % it.
        if  ~any(speciesName_list_idx_list(:) == speciesName_idx)
            speciesName_list_idx_list(length(speciesName_list_idx_list) + 1) = speciesName_idx;
        end
    end
    
    % We need to put the U back before N to match the other N order.
    N_idx = find(strcmp([speciesName_list{:}], "N"));
    U_idx = find(strcmp([speciesName_list{:}], "U"));
    speciesName_list_idx_list_idx = find(speciesName_list_idx_list==N_idx);
    speciesName_list_idx_list = horzcat([speciesName_list_idx_list(1:speciesName_list_idx_list_idx-1), ...
             U_idx,  ...
             speciesName_list_idx_list(speciesName_list_idx_list_idx:end)]);
    
    % speciesName_list_idx_list
    % count_list
    % popSSEM = cell2mat(count_list(1)) % original - wrong
    popSSEM = cell2mat(count_list(speciesName_list_idx_list(1)));
    for i = 2:length(speciesName_list_idx_list)
        idx = speciesName_list_idx_list(i);
        popSSEM = horzcat(popSSEM, cell2mat(count_list(idx)));
    end

end