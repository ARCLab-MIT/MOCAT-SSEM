function [scen_properties] = make_collision_pairs_MOCAT4(scen_properties, Su_species, S_species, D_species, N_species)
% This function is specific to MOCAT 4. The generalized function for other
% MOCAT versions is make_collision_pairs_SBM()
% This function takes a scen_properties object with a list of species and 
% the same species organized in species_cell into archetypical categories.
% It calculates and creates a set of species_pair objects which is stored
% at scen_properties.species_pairs and used to compile collision equations
% during model build.  The model is aware of trackability and maneuverability
% Object fragmentation counts are based on the NASA standard
% breakup model. 

%% Collision stuff
% Create matrix with pairs of species for collision modelling.
species_cross_pairs = nchoosek(scen_properties.species,2);
species_self_pairs = nchoosek(scen_properties.species,1);
species_self_pairs(:,2) = species_self_pairs(:,1);
species_pairs = vertcat(species_cross_pairs, species_self_pairs);

% Create an empty array of class species_pair_class
species_pairs_classes = species_pair_class.empty;

% Create array for fragments (used for n_f approach)
n_f = sym('n_f', [scen_properties.N_shell, 1]);

% For each species pair:
%   calculate and populate gammas: what is lost or gained from the species affected by the collision. We calculate these based on the tables in the introduction.
%   populate source_sinks : array of species affected by the collision.
% Save collision information in scen_properties.species_pairs: an array of species_pair_class instances. This class calculate phi and add the collision coefficients to the source sink equations.

% b = ProgressBar(numIterations,'Title', 'Progress');
for i = 1:height(species_pairs)
    % b(1, [], []);
    n1 = species_pairs(i,1).species_properties.sym_name;
    n2 = species_pairs(i,2).species_properties.sym_name;
    s1 = species_pairs(i,1);
    s2 = species_pairs(i,2);

    gammas = sym.empty;
    source_sinks = species.empty;

    % populate gammas and source_sinks for collisions that affect the
    % number of unslotted satellites
    if n1 == "Su" || n2 == "Su"
        gammas = horzcat(gammas, ones([scen_properties.N_shell, 1], 'sym'));
        index = width(gammas);
        source_sinks(index) = Su_species;
        
        if s1.species_properties.maneuverable && s2.species_properties.maneuverable
            gammas(:,index) = -s1.species_properties.alpha_active;
            if n1 == "S"
                gammas(:,index) = gammas(:,index) .* min(s2.species_properties.sym ./ (s1.species_properties.sym + s2.species_properties.sym),0);
            elseif n2 == "S"
                gammas(:,index) = gammas(:,index) .* min(s1.species_properties.sym ./ (s1.species_properties.sym + s2.species_properties.sym),0);
            end
        elseif n1 == "D" || n1 == "N" 
            gammas(:,index) = -(scen_properties.delta + s2.species_properties.alpha);
        elseif n2 == "D" || n2 == "N" 
            gammas(:,index) = -(scen_properties.delta + s1.species_properties.alpha);
        end
    end

    % populate gammas and source_sinks for collisions that affect the
    % number of unslotted satellites
    if n1 == "S" || n2 == "S"
        gammas = horzcat(gammas, ones([scen_properties.N_shell, 1], 'sym'));
        index = width(gammas);
        source_sinks(index) = S_species;
        
        if s1.species_properties.slotted && s2.species_properties.slotted
            gammas(:,index) = -s1.species_properties.alpha_active * (1-s1.species_properties.slotting_effectiveness);
        elseif n1 == "Su"
            gammas(:,index) = -s1.species_properties.alpha_active * min(s2.species_properties.sym ./ (s1.species_properties.sym + s2.species_properties.sym),0);
        elseif n2 == "Su"
            gammas(:,index) = -s2.species_properties.alpha_active * min(s1.species_properties.sym ./ (s1.species_properties.sym + s2.species_properties.sym),0);
        elseif n1 == "D" || n1 == "N" 
            gammas(:,index) = -(scen_properties.delta + s2.species_properties.alpha);
        elseif n2 == "D" || n2 == "N" 
            gammas(:,index) = -(scen_properties.delta + s1.species_properties.alpha);
        end
    end

    % populate gammas and source_sinks for collisions that affect the
    % number of derelict
    if (~s1.species_properties.maneuverable || ~s2.species_properties.maneuverable) && (n1 ~= "N" || n2 ~= "N")
        gammas = horzcat(gammas, ones([scen_properties.N_shell, 1], 'sym'));
        index = width(gammas);
        source_sinks(index) = D_species;
        if s1.species_properties.maneuverable
            gammas(:,index) = scen_properties.delta;
        elseif s2.species_properties.maneuverable
            gammas(:,index) = scen_properties.delta;
        else
            gammas(:,index) = -1;
        end
    end
    
    
    % Populate gammas and source_sinks for the debris species.
    index = width(gammas) + 1;
    gammas(:,index) = n_f;
    source_sinks(index) = N_species;

    if s1.species_properties.slotted && s2.species_properties.slotted
        gammas(:,index) = gammas(:,index) * s1.species_properties.alpha_active * (1-s1.species_properties.slotting_effectiveness);
    elseif s1.species_properties.maneuverable && s2.species_properties.maneuverable
        gammas(:,index) = gammas(:,index) * s1.species_properties.alpha_active;
    elseif s1.species_properties.maneuverable
        gammas(:,index) = gammas(:,index) * s1.species_properties.alpha;
    elseif s2.species_properties.maneuverable
        gammas(:,index) = gammas(:,index) * s2.species_properties.alpha;
    end

    % Make the class with the
    % equations.\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    species_pairs_classes(i) = species_pair_class(s1, s2, gammas, source_sinks, scen_properties);

    % Correct equations such that collisions between satellites and
    % derelict are catastrophic instead of damaging
    if (n1 == "S" || n1 == "Su" || n1 == "D") && (n2 == "S" || n2 == "Su" || n2 == "D")
        M1 = s1.species_properties.mass;
        M2 = s2.species_properties.mass;
        n_f_catastrophic = 0.1*scen_properties.LC^(-1.71)*(M1+M2)^(0.75)* ones(size(scen_properties.v_imp2));
        n_f_damaging = 0.1*scen_properties.LC^(-1.71)*(min(M1,M2)*scen_properties.v_imp2.^2).^(0.75); 
        for testi = 1:length(scen_properties.species) % find debris equation index
            if scen_properties.species(testi).species_properties.sym_name == "N"
                eq_index_N = testi;
                break
            end
        end
        N_eqn = species_pairs_classes(i).eqs(:,eq_index_N);
        species_pairs_classes(i).eqs(:,eq_index_N) = N_eqn.*transpose(n_f_catastrophic)./transpose(n_f_damaging);
    end
    
    % To print out species_pairs_classes for debugging (uncomment following lines)
    % species_pairs_classes(i)
    % species_pairs_classes(i).nf
    % species_pairs_classes(i).phi
    % gammas
    % T = vpa(species_pairs_classes(i).eqs,4)

end % End pairs loop
% b.release();
scen_properties.species_pairs = species_pairs_classes;
disp('Done')