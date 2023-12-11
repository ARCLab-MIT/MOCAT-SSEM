function [J,results] = capacity_objective(simulation, objective_settings)
    
    %% Fill in values for objective_settings

%     % Intrinsic Capacity
%     % sep_dist_method = 'distance' or 'angle'
%     sep_dist_method = 'distance';
%     sep_dist = 5; % minimum separation distance [km]
% 
%     %sep_dist_method = 'angle';
%     %sep_angle = 0.1; % minimum angular separation distance [deg]
%     
%     inc = 45; %degrees
% 
%     VAR.beta1 = 1e6; % Weight on the number of satellites.
%     VAR.beta3 = 1e1; % Weight on slotted vs. unslotted ratio constraint
% 
%     VAR.tol_debris = 30; % 0.1% of growth of debris population with respect to the initial conditions 
%     % VAR.tol_debris = max([sum(Ni)+1e4,5e4]); % absolute number of debris 
%     VAR.tol_debris_der = 20; % 0, tolerance of the derivative of the debris at the final time
%     VAR.threshold = 0.1;% threshold between slotted and unslotted satellites (not to cause too much difference)
%     VAR.S_Su = 0.1;% threshold between slotted and unslotted satellites (not to cause too much difference)
%     
    
    %% load results from intrinsic capacity and inclination vs altitude
    
    %TODO: this computes about a 0.2 second per-run penalty, but the
    %inclination really is an objective setting. Optimization could
    %computer this once and cache.
    
    filename = 'LawDF_min_500_fit_range_max_10000_sols_1.csv';
    A = readmatrix(filename);
    A = A(2:end,:);
    index_mat = A(:,1);
    inc_mat   = A(:,2); % [deg]
    c_mat     = A(:,3);
    b_mat     = A(:,4);
    R2_mat    = A(:,5);
    
    % Figure out which row to use based on the chosen inclination.
    [minValue,closestIndex] = min(abs(inc_mat-objective_settings.inc));
    ind_intrinsic = closestIndex; % it corresponds to 45 degree of inclination
 
    R02 = simulation.scen_properties.R02;
    re = simulation.scen_properties.re;
    rad = pi/180;
    if strcmp(objective_settings.sep_dist_method,'angle')
        N_sat_eq_ang = @(c,b) (objective_settings.sep_angle./c).^(1./b); % results from intrinsic capacity analysis
        N_sat = N_sat_eq_ang(c_mat(ind_intrinsic),b_mat(ind_intrinsic)).*(R02(2)-R02(1))./(sep_angle*rad*(R02(2:end)+re)); % N_sat_intrinsic*(n° shells per bin)
    else
        N_sat_eq = @(h,c,b) (objective_settings.sep_dist./(re+h)./rad./c).^(1./b); % results from intrinsic capacity analysis
        N_sat = N_sat_eq(R02(2:end),c_mat(ind_intrinsic),b_mat(ind_intrinsic))*(R02(2)-R02(1))/objective_settings.sep_dist; % N_sat_intrinsic*(n°  shells per bin)
    end
    N_sat_tot_intrinsic = sum(N_sat);

    
    %% propagation time constraint (if the final desired time is not reached)
    % TODO: Fix later.
%     if T(end)<VAR.tf % propagation doesn't reach the desired final time
%         pen_time = 1e9*(VAR.tf-T(end));
%     else
%         pen_time = 0;
%     end
    pen_time = 0;
    %% Read results

    S_list = get_species_name_list_from_species_cell(simulation, "S");
    Su_list = get_species_name_list_from_species_cell(simulation, "Su");
    UN_list = [get_species_name_list_from_species_cell(simulation, "U")  get_species_name_list_from_species_cell(simulation, "N")];

    S_results = get_species_results_from_species_names(simulation, S_list);
    S_results = sum(S_results, 3); % In case more than one species
    Su_results = get_species_results_from_species_names(simulation, Su_list);
    Su_results = sum(Su_results, 3); % In case more than one species
    UN_results = get_species_results_from_species_names(simulation, UN_list);
    UN_results = sum(UN_results, 3); % In case more than one species
    
    N_tot_real_S = sum(S_results(end,:)); % total number of slotted satellites 
    N_tot_real_Su = sum(Su_results(end,:)); % total number of unslotted satellites 
    N_tot_real = N_tot_real_S+N_tot_real_Su; % total number of satellites
    N_der = diff(UN_results(end-1:end,:))/(simulation.results.T(end)-simulation.results.T(end-1));
    %% intrinsic capacity constraint

    diff_S = N_sat-S_results;
    pen_sat = 0;
    pen_const = 0;
    for i=1:size(S_results,1) % check for intrinsic capacity constraint for each altitude and for each times step
        pen_sat_it_neg = diff_S(i,diff_S(i,:)<0);
        if ~isempty(pen_sat_it_neg)
            pen_sat = pen_sat + sum(abs(pen_sat_it_neg)); % it penalizes when intrinsic capacity is exceeded
            pen_const = pen_const + length(simulation.scen_properties.scen_times)/i; % it penalizes the time step when intrinsic capacity is exceeded
        end
    end
    %% final debris population constraint
    
    N_der_pos = N_der>objective_settings.tol_debris_der;
    if isempty(N_der_pos)
        pen_debris = 0;
    else
        pen_debris = 100*sum((N_der(N_der_pos))); % 100
    end
    
    %% slotted vs unslotted ratio constraint
    
    check_S_Su = Su_results(end,:)-objective_settings.threshold*S_results(end,:);
    check_S_Su_excess = check_S_Su(check_S_Su>0);
    
    if ~isempty(check_S_Su_excess) 
        pen_unsl = objective_settings.beta3*sum(check_S_Su_excess);
    else
        pen_unsl = 0;
    end
    
    %% Compile Objective.
    J = objective_settings.beta1/N_tot_real+pen_sat+pen_const+pen_time+pen_debris+pen_unsl;
    
    results{1} = simulation.results.T;
    results{2} = simulation.results.X;
    results{3} = [pen_time pen_unsl pen_sat+pen_const pen_debris];
end
