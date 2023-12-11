function [scen_properties] = make_collision_pairs_SBM(scen_properties, N_species)
% This function takes a scen_properties object with a list of species and 
% the same species organized in species_cell into archetypical categories.
% It calculates and creates a set of species_pair objects which is stored
% at scen_properties.species_pairs and used to compile collision equations
% during model build.  The model is aware of trackability, maneuverability,
% and slotting. Object fragmentation counts are based on the NASA standard
% breakup model.

%% Collision stuff
species_cross_pairs = nchoosek(scen_properties.species,2);

species_self_pairs = nchoosek(scen_properties.species,1);
species_self_pairs(:,2) = species_self_pairs(:,1);
species_pairs = vertcat(species_cross_pairs, species_self_pairs);
species_pairs_classes = species_pair_class.empty;
n_f = sym('n_f', [scen_properties.N_shell, 1]);

if nargin < 2
    N_species = scen_properties.species_cell.N;
end

% Compute the Mass bin centers, edges, and widths once rather for all
% debris species for fragmentation generation code.
binC =     zeros(1, length(N_species)); %Bin centers
binE =     zeros(1, 2*length(N_species)); % Bin edges (unused
binW =     zeros(1, length(N_species)); % Bin widths
LBgiven = scen_properties.LC;

for N_index = 1:length(N_species)
    % disp([num2str(N_index), " ", num2str(2*N_index-1), " ", num2str(2*N_index)])
    binC(N_index) = N_species(N_index).species_properties.mass;
    binE(2*N_index-1:2*N_index) = [N_species(N_index).species_properties.mass_lb, N_species(N_index).species_properties.mass_ub];
    binW(N_index) = N_species(N_index).species_properties.mass_ub - N_species(N_index).species_properties.mass_lb;
end

binE = unique(binE); % Combine connected lower and upper bounds

numIterations = height(species_pairs);
disp("Creating species pairs...")
b = ProgressBar(numIterations,'Title', 'Progress');
for i = 1:height(species_pairs)
    b(1, [], []);
    n1 = species_pairs(i,1).species_properties.sym_name;
    n2 = species_pairs(i,2).species_properties.sym_name;
    s1 = species_pairs(i,1);
    s2 = species_pairs(i,2);
    gammas = -ones([scen_properties.N_shell, 2], 'sym'); % Seed with what is lost from each of the two colliding species
    source_sinks = [s1, s2];
    % b.printMessage(sprintf('@Iteration %i' + " Species1:" + n1 + " T:" + string(s1.species_properties.trackable) ...
    %                                                           + " M:" + string(s1.species_properties.maneuverable)...
    %                                                           + " S:" + string(s1.species_properties.slotting_effectiveness)...
    %                                        + "     Species2:" + n2 + " T:" + string(s2.species_properties.trackable) ...
    %                                                           + " M:" + string(s2.species_properties.maneuverable)...
    %                                                           + " S:" + string(s2.species_properties.slotting_effectiveness), i));
    % 
    % Figure out reductions to primary species from col.
    % alpha_a * alpha_a if both maneuverable
    if s1.species_properties.maneuverable && s2.species_properties.maneuverable
        gammas(:,1) = gammas(:,1) * s1.species_properties.alpha_active * s2.species_properties.alpha_active;
        % (1-zeta) using worst zeta if both are slotted
        if s1.species_properties.slotted && s2.species_properties.slotted
            gammas(:,1) = gammas(:,1) * min(s1.species_properties.slotting_effectiveness, s2.species_properties.slotting_effectiveness);
        end
    % alpha is one is maneuverable and both are trackable.
    elseif (s1.species_properties.maneuverable && ~s2.species_properties.maneuverable) || (s2.species_properties.maneuverable && ~s1.species_properties.maneuverable)
        %TODO: If a species is maneuverable, but not trackable, this will
        % not produce a reduction.
        if s1.species_properties.trackable && s2.species_properties.maneuverable
            gammas(:,1) = gammas(:,1) * s2.species_properties.alpha;
        elseif s2.species_properties.trackable && s1.species_properties.maneuverable
            gammas(:,1) = gammas(:,1) * s1.species_properties.alpha;   
        end
    end
    % The gamma burden is symmetric lost to both colliding species.
    gammas(:,2) = gammas(:,1);
    
    %     % Add the fragments that get made
    %     for index = 1:length(N_species)
    %         gammas(:,2+index) = 1 * n_f;
    %         source_sinks(2+index) = N_species(index);
    %     end %End debris type loop

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
    end
    
    %disp("    Gamma1: " + string(gammas(1,1)))

    % Populate gammas and source_sinks for the debris species.
    for index = 1:length(N_species)
        gammas(:,2+index) = -gammas(:,1) .* frags_made(:, index); % Because frags only happen in collisions.
        source_sinks(2+index) = N_species(index);
    end %End debris type loop

    % Make the class with the equations.

    species_pairs_classes(i) = species_pair_class(s1, s2, gammas, source_sinks, scen_properties);
    
    % To see resultant equations associated with species_pairs_classes(i)
    %for spec_index = 1:length(species_list)
    %    spec_names(spec_index) = species_list(spec_index).species_properties.sym_name;
    %end
    %T = array2table(species_pairs_classes(i).eqs, 'VariableNames', spec_names);

    
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
scen_properties.species_pairs = species_pairs_classes;
disp("Done")