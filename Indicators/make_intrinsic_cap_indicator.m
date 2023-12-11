
function ind_struct = make_intrinsic_cap_indicator(scen_properties, sep_dist_method, varargin)
    
    % Intrinsic Capacity
    % sep_dist_method = 'distance' or 'angle'
    % Default sep_dist is 25 km for distance 
    % Default sep_angle is 0.2 degree for angle
    % default_shell_sep is 5 km
    % Angle used for intrinsic capacity defaults to 45 degrees (optimistic case)    
    
    %TODO: Default to
    %https://www.space-track.org/documents/Spaceflight_Safety_Handbook_for_Operators.pdf
    % Table 3.
    
    default_inc = 45.0; %degrees [Supported range is 0<x<=90. For retrograde orbits, use prograde equivalent (e.g. 94 deg = 86 deg)
    default_sep_angle = 0.2; % minimum in-shell angular separation distance [deg]
    default_sep_distance = 25.0; % minimum in-shell angular separation distance [km]
    default_shell_sep = 5; % average assumed exclusive height of orbital shells. [km]
    default_graph = false; % bool for whether to plot intrinsic capacity per bin

    inc_validation = @(x) isnumeric(x) && isscalar(x) && x > 0 && x<=90 ;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && x>0;
    p = inputParser;
    addRequired(p, 'scen_properties', @(x) class(x) == "scen_properties_class")
    addRequired(p, 'sep_dist_method', @(x) any(strcmp(x, ["angle", "distance"])))
    addOptional(p,'sep_angle',default_sep_angle, validScalarPosNum);
    addOptional(p,'sep_dist',default_sep_distance, validScalarPosNum);
    addOptional(p,'shell_sep',default_sep_distance, validScalarPosNum);
    addOptional(p, 'inc', default_inc, inc_validation);
    addOptional(p, 'graph', default_graph, @(x) islogical(x));

    parse(p, scen_properties, sep_dist_method, varargin{:})
    % Provide warning for mismatched sep_dist_method and distance value.
    if strcmp(sep_dist_method, "angle") && ~ismember(p.UsingDefaults, 'sep_dist')
        warning("'sep_dist' argument was passed despite 'sep_dist_method' " + ...
                "being set to 'angle'. 'sep_dist' parameter was ignored. " + ...
                "'sep_angle' may have been set to a default value") 
    end

    if strcmp(sep_dist_method, "distance") && ~ismember(p.UsingDefaults, 'sep_angle')
        warning("'sep_angle' argument was passed despite 'sep_dist_method' " + ...
                "being set to 'distance'. 'sep_angle' parameter was ignored. " + ...
                "'sep_dist' may have been set to a default value") 
    end

    shell_sep = p.Results.shell_sep;
    inc = p.Results.inc;
    if isfield(p.Results, 'sep_angle')
        sep_angle = p.Results.sep_angle;
    end

    if isfield(p.Results, 'sep_dist')
        sep_dist = p.Results.sep_dist;
    end

    filename = 'LawDF_min_500_fit_range_max_10000_sols_1.csv';
    A = readmatrix(filename);
    A = A(2:end,:);
    index_mat = A(:,1);
    inc_mat   = A(:,2); % [deg]
    c_mat     = A(:,3);
    b_mat     = A(:,4);
    R2_mat    = A(:,5);
    
    % Figure out which row to use based on the chosen inclination.
    [~,closestIndex] = min(abs(inc_mat-inc));
    ind_intrinsic = closestIndex;
    
    R02 = scen_properties.R02;
    re = scen_properties.re;
    rad = pi/180;
    if strcmp(sep_dist_method,'angle')
        N_sat_eq_ang = @(c,b) (sep_angle./c).^(1./b); % results from intrinsic capacity analysis
        % Original approach
        % N_sat = N_sat_eq_ang(c_mat(ind_intrinsic),b_mat(ind_intrinsic)).*(R02(2)-R02(1))./(sep_angle*rad*(R02(2:end)+re)); % N_sat_intrinsic*(n° shells per bin)
        shells_per_bin = diff(R02(1:end))./shell_sep;
        N_sat = N_sat_eq_ang(c_mat(ind_intrinsic),b_mat(ind_intrinsic)).* shells_per_bin ; % N_sat_intrinsic*(n° shells per bin)
    else
        N_sat_eq = @(h,c,b) (sep_dist./(re+h)./rad./c).^(1./b); % results from intrinsic capacity analysis
        shells_per_bin = diff(R02(1:end))./shell_sep;
        N_sat = N_sat_eq(R02(2:end),c_mat(ind_intrinsic),b_mat(ind_intrinsic)).*shells_per_bin; % N_sat_intrinsic*(n°  shells per bin)
    end

    % N_sat_tot_intrinsic = sum(N_sat);
    
    slotted_sats_eqs = sym(zeros(scen_properties.N_shell,1));
    for species_i = 1:length(scen_properties.species)
        species = scen_properties.species(species_i);
        if species.species_properties.slotted == true
            % s_name = species.species_properties.sym_name;
            slotted_sats_eqs = slotted_sats_eqs +species.species_properties.sym;
        end
    end
    
    unconsumed_intrinsic_capacity = N_sat.'-slotted_sats_eqs;
    ind_struct = indicator_var_class("unconsumed_intrinsic_capacity", "manual", [], unconsumed_intrinsic_capacity);

    if p.Results.graph == true
        figure()
        h = scatter(R02(2:end).', N_sat.');
        title("Intrinsic Capacity per Altitude Bin")
        xlabel("Altitude of Bin Ceiling [km]") 
        ylabel("Intrinsic Capacity [Satellites]") 
        disp("Done")
    end
