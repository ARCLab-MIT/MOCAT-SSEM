function [pairs, species_names, table] = generate_species_pairs(scen_properties, S_species, D_species, N_species)
    
    species_cross_pairs = nchoosek(scen_properties.species,2);
    species_self_pairs = nchoosek(scen_properties.species,1);
    species_self_pairs(:,2) = species_self_pairs(:,1);
    species_pairs = vertcat(species_cross_pairs, species_self_pairs);
    species_pairs_classes = species_pair_class.empty;
    
    n_f = sym('n_f', [scen_properties.N_shell, 1]);
    
    % Compute the Mass bin centers, edges, and widths once rather for all
    % debris species for fragmentation generation code.
    binC =     zeros(1, length(N_species)); %Bin centers
    binE =     zeros(1, 2*length(N_species)); % Bin edges (unused
    binW =     zeros(1, length(N_species)); % Bin widths
    LBgiven = scen_properties.LC;
    
    for N_index = 1:length(N_species)
        %disp([num2str(N_index), " ", num2str(2*N_index-1), " ", num2str(2*N_index)])
        binC(N_index) = N_species(N_index).species_properties.mass;
        binE(2*N_index-1:2*N_index) = [N_species(N_index).species_properties.mass_lb, N_species(N_index).species_properties.mass_ub];
        binW(N_index) = N_species(N_index).species_properties.mass_ub - N_species(N_index).species_properties.mass_lb;
    end
    
    binE = unique(binE); % Combine connected lower and upper bounds
    
    numIterations = height(species_pairs);
    % Create species pairs
    b = ProgressBar(numIterations,'Title', 'Progress');
    for i = 1:height(species_pairs)
        b(1, [], []);
        n1 = species_pairs(i,1).species_properties.sym_name;
        n2 = species_pairs(i,2).species_properties.sym_name;
        s1 = species_pairs(i,1);
        s2 = species_pairs(i,2);
    
        gammas = [];
        source_sinks = species.empty;
        % populate gammas and source_sinks for collisions that affect the
        % number of unslotted satellites
        if n1 == "S" || n2 == "S"
            gammas = horzcat(gammas, ones([scen_properties.N_shell, 1], 'sym'));
            index = width(gammas);
            source_sinks(index) = S_species;
            
            if n1 == "S" && n2 == "S"
                gammas(:,index) = -s1.species_properties.alpha_active;
            elseif n1 == "D" || n1 == "N" 
                gammas(:,index) = -(scen_properties.delta + s2.species_properties.alpha_active);
            elseif n2 == "D" || n2 == "N" 
                gammas(:,index) = -(scen_properties.delta + s1.species_properties.alpha_active);
            end
        end
    
        % populate gammas and source_sinks for collisions that affect the
        % number of derelict
        if (n1 ~= "S" || n2 ~= "S") && (n1 ~= "N" && n2 ~= "N")
            gammas = horzcat(gammas, ones([scen_properties.N_shell, 1], 'sym'));
            % source_sinks = horzcat(horzcat, D_species);
            index = width(gammas);
            source_sinks(index) = D_species;
            if n1 == "S"
                gammas(:,index) = scen_properties.delta;
            elseif n2 == "S"
                gammas(:,index) = scen_properties.delta;
            else
                gammas(:,index) = -1;
            end
        end
    
        % Find debris generation for each debris class from S1-S2 collision.
        m1 = s1.species_properties.mass;
        m2 = s2.species_properties.mass;
        r1 = s1.species_properties.radius;
        r2 = s2.species_properties.radius;
        RBflag = max(s1.species_properties.RBflag, s2.species_properties.RBflag); % True if either is rocket body
        
        % Initialize empty array. Rows = altitudes with different dv values,
        % cols = debris species in order of debris species list..
        frags_made = zeros(length(scen_properties.v_imp2), length(N_species));
        isCatastrophic = zeros(1, length(N_species));
    
        % Populate frags_made based on created fragments for each bin and
        % debris mass.
        for dv_index = 1:length(scen_properties.v_imp2)
            dv = scen_properties.v_imp2(dv_index);
            %[nums, isCatastrophic, binOut]  = EVOLVEbins(m1,m2,r1,r2,dv,[],binE,[], LBgiven, RBflag)
            [frags_made(dv_index, :), isCatastrophic(1,dv_index), ~] = EVOLVEbins(m1,m2,r1,r2,dv,[],binE,[], LBgiven, RBflag);
            % The number of fragments is computed by the EVOLVEbins function
        end
        
        start_index = width(gammas);
        
        % Populate gammas and source_sinks for the debris species.
        for index = 1:length(N_species)
            gammas(:,start_index+index) = frags_made(:, index);
            source_sinks(start_index+index) = N_species(index);
    
            if s1.species_properties.maneuverable && s2.species_properties.maneuverable
                gammas(:,start_index+index) = gammas(:,start_index+index) * s1.species_properties.alpha_active;
            elseif s1.species_properties.maneuverable
                gammas(:,start_index+index) = gammas(:,start_index+index) * s1.species_properties.alpha;
            elseif s2.species_properties.maneuverable
                gammas(:,start_index+index) = gammas(:,start_index+index) * s2.species_properties.alpha;
            end
        end %End debris type loop
            
    
        % Make the class with the equations.
    
        species_pairs_classes(i) = species_pair_class(s1, s2, gammas, source_sinks, scen_properties);
        clear spec_names
        % To see resultant equations associated with species_pairs_classes(i)
        species_list = scen_properties.species;
        for spec_index = 1:length(species_list)
           spec_names(spec_index) = species_list(spec_index).species_properties.sym_name;
        end
        T = array2table(species_pairs_classes(i).eqs, 'VariableNames', spec_names);
        disp(T); % Labeled table of initial populations of each species
    
        
        % For equation comp testing:
        %delta = sym('delta');
        %zeta = sym('zeta'); 
        %alpha = sym('alpha');
        %alpha_a = sym('alpha_a');
        %     for source_sink_i =1:width(source_sinks)
        %         source_sinks_names(source_sink_i) = source_sinks(source_sink_i).species_properties.sym_name; 
        %    end
    end % End pairs loop

b.release();
pairs = species_pairs_classes;
species_names = spec_names;
table = T;