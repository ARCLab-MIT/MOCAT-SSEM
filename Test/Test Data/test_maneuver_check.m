%% Make Scenario Properties
clear
clc

fileName = matlab.desktop.editor.getActiveFilename;
folder = fileparts(which(fileName)); 
addpath(genpath(folder));
cd(folder);
addpath(genpath(pwd))
addpath(genpath(folder));
VAR = MOCATSSEM_VAR_Cons();


% Add the time stuff needed
VAR.start_year = 2022;
VAR.start_month = 12;
VAR.start_day = 1;
VAR.start_hour = 0;
VAR.start_minute = 0;
VAR.start_second = 0;
VAR.start_second_fraction = 0;
VAR.simulation_duration = 200.0; %Years
VAR.N_steps = 200;
VAR.scen_times = linspace(0,VAR.simulation_duration,VAR.N_steps)';

VAR.density_filepath = fullfile(".", "Atmosphere Model", "JB2008", "Precomputed", "dens_highvar_2000.mat");
VAR.dens_model = @JB2008_dens_func;
VAR.time_dep_density = true;

%VAR.dens_model = @static_exp_dens_func;
%VAR.time_dep_density = false;

VAR.integrator = @ode15s; % Fixed options: ode1, ode2, ode3, ode4, ode5
%                                     % Variable options: ode45, ode21,
%                                     % ode113, oder78, ode89, ode15s,
%                                     % oder23s, ode23t, ode23tb
VAR.integrator_type = 'variable'; %or 'fixed';

scen_properties = scen_properties_class(VAR);

%% Species Properties
species_properties = {};

% Capabilities
species_properties.drag_effected = false; % bool, does object decrease altitude due to drag
species_properties.active = true; % bool, use alpha_active for collisions vs. trackable
species_properties.maneuverable = true; % bool, use alpha_active for collisions vs. trackable
species_properties.trackable = true; % bool, others can avoid this with alpha
species_properties.RBflag = 0; % No Rocket Bodies yet.


%% Make Species
t = 10;
my_drag_func = @drag_func;
null_launch_func = @launch_func_null;
launch_func = @launch_func_gauss;
pmd_func = 0;

% S
species_properties.sym_name = "S";
species_properties.Cd = 2.2; % unitless
species_properties.mass = [255]; % kg
species_properties.radius = [sqrt(15.45/pi)]; % m
species_properties.A = species_properties.radius.^2 * pi; % m^2 %TODO don't guess area.
% For this amr and others, we use first mass so multi_property_species will
% scale it automatically for us using the scalred radii.
species_properties.amr = species_properties.A./species_properties.mass;  % m^2/kg
species_properties.beta = species_properties.Cd*species_properties.amr; % ballistic coefficient

% Orbit Properties
species_properties.slotted = false; % bool
species_properties.slotting_effectiveness  = 1.0; % Float [0,1],

% Capabilities
species_properties.drag_effected = false; % bool, does object decrease altitude due to drag
species_properties.active = true; % bool, use alpha_active for collisions vs. trackable
species_properties.maneuverable = true; % bool, use alpha_active for collisions vs. trackable
species_properties.trackable = true; % bool, others can avoid this with alpha.
species_properties.deltat = [6]; % lifetime of spacecraft in years
species_properties.Pm = .90; % Post-mission disposal efficacy Float [0,1]
species_properties.alpha = 0; % Efficacy of Collision Avoidance vs. inactive
species_properties.alpha_active = 0; % Efficacy of Collision Avoidance vs. other active 
species_properties.orbit_raising = 'onstation'; % Whether orbit-raising will occur.
species_properties.insert_altitude = 250.0; % Altitude in km for insertion altitude
species_properties.onstation_altitude = 7500.0; % Altitude in km for onstation altitude
species_properties.RBflag = 0; % No Rocket Bodies yet.
species_properties.lambda_constant = 2000;
S_species = species(@launch_func_constant, @pmd_func_sat, my_drag_func, species_properties, scen_properties);

% Sns 3U cubesat and 250 kg sat.
species_properties.sym_name = "Sns";
species_properties.Cd = 2.2; % unitless
species_properties.deltat = [3 ];
species_properties.mass = [255]; % kg
species_properties.radius = [sqrt(15.45/pi) ]; % m
species_properties.A = species_properties.radius.^2 * pi; % m^2 %TODO don't guess area.
species_properties.amr = species_properties.A./species_properties.mass;  % m^2/kg
species_properties.beta = species_properties.Cd*species_properties.amr; % ballistic coefficient
species_properties.slotted = false; % bool
%species_properties.lambda_constant = 5000;
Sns_species = species(@launch_func_constant, @pmd_func_sat, my_drag_func, species_properties, scen_properties);


% N
% Assuming al spheres of 1 cm, 10cm diameter, 15 kg 
species_properties.sym_name = "N";
species_properties.Cd = 2.2; % unitless
species_properties.deltat = NaN;
species_properties.mass =   [0.00141372     0.5670     255]; 
species_properties.radius = [0.01           0.1321    sqrt(15.45/pi)]; % m
species_properties.A = pi*species_properties.radius.^2;
species_properties.amr = species_properties.A./species_properties.mass;  % m^2/kg
species_properties.beta = species_properties.Cd*species_properties.amr; % ballistic coefficient
species_properties.drag_effected = true; % bool, does object decrease altitude due to drag
species_properties.maneuverable = false; % bool, use alpha_active for collisions vs. trackable
species_properties.active = false;
N_species = multi_property_species(null_launch_func, @pmd_func_none,my_drag_func, species_properties, scen_properties).species_list;

% % N
% % Assuming al spheres of 1 cm, 10cm diameter, 15 kg 
% species_properties.sym_name = "N";
% species_properties.Cd = 2.2; % unitless
% species_properties.deltat = NaN;
% species_properties.mass =   [0.1  1.5   255]; 
% species_properties.radius = [0.01 0.5   sqrt(15.45/pi) ]; % m
% species_properties.A = pi*species_properties.radius.^2;
% species_properties.amr = species_properties.A./species_properties.mass;  % m^2/kg
% species_properties.beta = species_properties.Cd*species_properties.amr; % ballistic coefficient
% species_properties.drag_effected = true; % bool, does object decrease altitude due to drag
% species_properties.maneuverable = false; % bool, use alpha_active for collisions vs. trackable
% species_properties.active = false;
% N_species = multi_property_species(null_launch_func, @pmd_func_none,my_drag_func, species_properties, scen_properties).species_list;

%Graft on the debris properties.
start_i = 1+ length(species_properties.mass) - length(S_species) -  length(Sns_species);
end_i = length(N_species);
linked_spec_list = [S_species Sns_species];

% Assign matching debris increase for a species due to failed PMD
for i = 1:length(linked_spec_list)
    found_mass_match_debris = false; 
    spec_mass = linked_spec_list(i).species_properties.mass;
    %disp(linked_spec_list(i).species_properties.sym_name + " " + string(spec_mass))
    for species_i = 1:length(N_species)
        if spec_mass == N_species(species_i).species_properties.mass
            %disp("   N_mass " + string(N_species(species_i).species_properties.mass))
            N_species(species_i).pmd_func = @pmd_func_derelict;
            N_species(species_i).species_properties.pmd_linked_species = [N_species(species_i).species_properties.pmd_linked_species; linked_spec_list(i)];
            disp("Matched species " ...
                + linked_spec_list(i).species_properties.sym_name ...
                + " to debris species " ...
                +  N_species(species_i).species_properties.sym_name + ".")
            found_mass_match_debris = true;
        end
    end
    if ~found_mass_match_debris
        warning("No matching mass debris species found for species " + ...
                linked_spec_list(i).species_properties.sym_name + ...
                " with mass " + string(spec_mass) + "." )
    end
end

species_list = [S_species, N_species, Sns_species];
scen_properties.species = species_list;
scen_properties.species_cell = struct('S', S_species, 'N', N_species, 'Sns', Sns_species);

%% Starting Pop and FLM
% Display starting population
for spec_index = 1:length(scen_properties.species)
   spec_names(spec_index) = scen_properties.species(spec_index).species_properties.sym_name;
end

%load("ws_ADEPT_Testing_Take2.mat")
[x0,FLM_steps] = ADEPT_Traffic_Model("D:\ADEPT Data\start_full.asem.csv", scen_properties);

disp("Done loading initial population")

if length(unique(round(diff(scen_properties.scen_times),5))) == 1
    time_step = unique(round(diff(scen_properties.scen_times),5));
else
    error("VnV not stet up for variable time step runs.")
end

for ii = 1:length(scen_properties.species)
    sel_species = scen_properties.species(ii);
    species_names = sel_species.species_properties.sym_name;
    %spec_index = get_species_index_from_name(scen_properties, species_names);

    spec_FLM = array2table(zeros(scen_properties.N_shell,length(scen_properties.scen_times)), ...
                       RowNames = string(1:scen_properties.N_shell), ...
                       VariableNames = string(scen_properties.scen_times));
    
    % Make the FLM array (Nshell x NTimes)
    % Assume no launches at Tstep 1, as the pop is considered x0
    for jj = 2:length(FLM_steps) % JJ iterates over times
        spec_FLM(:,jj) = FLM_steps{jj}(:,ii); % ii is species
    end
    
    % Scale sats to sats/year and make an interp
    Lambdadot = cell(size(scen_properties.N_shell,2), 0);
    for hh = 1:scen_properties.N_shell
        x = scen_properties.scen_times;
        v = spec_FLM{hh,:}.'/time_step; % Convert from sats to sats/year
        Lambdadot{hh, 1} = @(X, t) interp1(x, v, t, 'linear');
    end
    sel_species.species_properties.lambda_funs = Lambdadot;
    sel_species.launch_func = @launch_func_lambda_fun;
end

disp("Done assigning future launch traffic model.")

%% Collision stuff
[scen_properties] = make_collision_pairs_SBM(scen_properties);
disp("Done")

%% Make AG SEM Goal Indicators

% LTS
indicator_var_list = indicator_var_class.empty();
% TODO: Hard because need both derivative and total.

% Short-Term Collision Avoidance
percentage = false;
per_species = false;
short_term_ind = make_active_loss_per_shell_counter(scen_properties, percentage, per_species);
all_col_indicators = make_collisions_per_species(scen_properties);

per_spacecraft = false;
ca_man_struct = make_ca_counter(scen_properties, "maneuverable", "trackable", "per_species", true, "per_spacecraft", per_spacecraft);

% Orbital Volume
intrinsic_ind = make_intrinsic_cap_indicator(scen_properties, "distance", ...
                                                    'sep_dist', 60, ...
                                                    'inc', 40, 'shell_sep', 5, ...
                                                    'graph', false);

indicator_var_list = [short_term_ind; all_col_indicators; ca_man_struct; intrinsic_ind; ];
scen_properties.indicator_var_list = indicator_var_list;
%%
x0_run = reshape(table2array(x0), [], 1);% Pad initial conditions with the indicator variables
my_sim = simulation_class(species_list, scen_properties);
my_sim.build_model();
my_sim.scen_properties.integrator = @ode15s;
my_sim.scen_properties.disp_times = false;


%%
my_sim.run_model(x0_run, 'progressBar', true, 'disp_times', false, 'simulation_duration',200, 'N_step', 2); %'start_date', today('datetime'), 
my_sim.pop_vs_time_vis()
disp("Done")
%%
% Intrinsic Capacity
my_sim.per_shell_indicator("unconsumed_intrinsic_capacity", "constraint", 0.0, 'constraint_type', "lower")

% Operational CA

%my_sim.per_shell_indicator("active_aggregate_collisions", "constraint", 1.0)

for ii = 1:length(scen_properties.indicator_var_list)
    disp(scen_properties.indicator_var_list(ii).name)
    if contains(scen_properties.indicator_var_list(ii).name, "aggregate_collisions")
        my_sim.per_shell_indicator(scen_properties.indicator_var_list(ii).name, "constraint", 1e-4, 'constraint_type', "upper")
    end
end

for ii = 1:length(scen_properties.indicator_var_list)
    disp(scen_properties.indicator_var_list(ii).name)
    if contains(scen_properties.indicator_var_list(ii).name, "maneuvers")
        my_sim.per_shell_indicator(scen_properties.indicator_var_list(ii).name, "constraint", 1, 'constraint_type', 'upper')
    end
end

% Long-Term Sustainability
my_sim.total_deb_species_deriv_evolv_vis("percentage", true, "color_scale_crop_percentiles", [5, 95], "constraint", 5)
my_sim.total_deb_pop_time_deriv("constraint", 5)

%% 
% Compare

% filepaths = dir("D:\ADEPT Data\results\scenario_case-A0\*\snapshot*species.csv");
% % iterate over folder
% 
% data = {};
% 
% % bin results for each snapshot into species, 
% for ii = 1:numel(filepaths)
%     disp(string(ii) + " / " + string(numel(filepaths)))
%     file_path = fullfile(filepaths(ii).folder, filepaths(ii).name);
% 
%     mc = split(filepaths(ii).name,"mc"); mc = mc{2}; mc = split(mc,"_year");
%     mc = mc{1};
% 
%     year = split(filepaths(ii).name,"_year"); year = year{2}; year = split(year,".txt");
%     year = str2num(year{1});
% 
%     data.("mc_" + mc).("year_" + year) =  bin_mass_shell_alt(file_path, scen_properties);
% end

load("ws_ADEPT_Testing_Take1.5 (mod2 results).mat", "data")

%%
% % Add average for each MC.
% 
% years_fields = fieldnames(data.mc_001);
% fn = fieldnames(data);
% for j = 1:numel(years_fields)
%     year = string(years_fields{j})
%     for k = 1:numel(fn)
%         field = fn{k}
%         if contains(field, "mc")
%             if k == 1
%                 data.(year)  = data.(field).(year);
%             else
%                 data.(year) = data.(year) + data.(field).(year);
%             end
%         else
%             disp("Skipping field " + field);
%         end
%     end
%     data.(year) = data.(year)./numel(fn);
% end
% 
% %load("ADEPT_Binned_result_data_by_year")
% disp("Ready to Compare")
% 
% my_sim.results.species_array

%%
match_years = [2024, 2025, 2026, 2027, 2028, 2029, 2030, 2031, 2032, 2033, 2034, 2035, 2043, 2045, 2050, 2055, 2065, 2075, 2100, 2125, 2150, 2200];
match_indexs = [];
for ii = 1:numel(match_years)
    match_year = match_years(ii);
    sim_dates = my_sim.scen_properties.start_date + years(my_sim.results.T);
    [min_dist, min_i] = min(abs(sim_dates - datetime(match_year, 1, 1)));
     match_indexs(ii) = min_i;
end
%%
% figure()
% for ii = 1:numel(match_indexs)
%     match_index = match_indexs(ii);
%     sim_results = squeeze(my_sim.results.species_array(match_index, : , :));
%     ADEPT_results = data.("year_" + string(match_years(ii)));
%     diff = table2array(ADEPT_results - sim_results);
%     alts = my_sim.scen_properties.HMid.';
%     % 
%     % surf(alts,categorical(spec_names).',diff.','edgecolor','none', 'LineStyle','none')
%     % %surf(alts,categorical(spec_names).',diff.','edgecolor','none', 'LineStyle','none')
%     % c = colorbar;
%     % caxis([-2500 2500]);
%     % c.Label.String = 'Difference in Object Count';
%     % zlabel('Difference in Object Count')
%     % xlabel('Altitude (km)')
%     % ylabel('Species')
%     % xlim([min(alts) max(alts)])
%     % zlim([-2500 2500])
%     % title('ADEPT Population (average) - MOCAT-SSEM for year = ' + string(match_years(ii)), 'Interpreter', 'none')
%     % grid on
%     % %view (0, 90)
% 
%     agg_diffs = sum(diff, 1);
%     % figure()
%     %c = colorbar;
%     %caxis([-2500 2500]);
%     %c.Label.String = 'Difference in Object Count';
%     scatter(categorical(spec_names), abs(agg_diffs))
%     xlabel('Species')
%     ylabel('Difference in Object Count')
%     title('ADEPT Population (average) - MOCAT-SSEM for year = ' + string(match_years(ii)), 'Interpreter', 'none')
%     grid on
%     set(gca, 'Yscale', 'log')
%     ylim([1e-2 1e10])
% end


%%
% Compare data by Species
alts = my_sim.scen_properties.HMid;

ADEPT_years = years(datetime(match_years, 12, 31) - my_sim.scen_properties.start_date);
%ADEPT_results = data.("year_" + string(match_years(ii))

for i = 1:length(my_sim.species_list)   
    figure;
    species_array = my_sim.results.species_array;
    surf(my_sim.results.T,alts,transpose(species_array(:,:,i)), 'edgecolor','none', "DisplayName", "MOCAT-SSEM", "FaceAlpha", 0.1, "EdgeAlpha", 0.1)
    hold on
    c = colorbar;
    c.Label.String = 'Number of Objects';
    zlabel('Number of Objects')
    xlabel('Time (years)')
    ylabel('Altitude (km)')
    xlim([0 max(my_sim.results.T)])
    ylim([min(alts) max(alts)])
    title('Population of Species '+strrep(my_sim.species_list(i).species_properties.sym_name,'p','.'), 'Interpreter', 'none')
    grid on

    ADEPT_results = zeros(numel(alts), numel(match_years)); % Alts by years for species
    for ii = 1:length(match_years)
        ADEPT_results_year = data.("year_" + string(match_years(ii))); % Alts x Species in fixed year
        ADEPT_results(:,ii) = ADEPT_results_year{:,i};
        if ii == 1
            scatter3( repelem(ADEPT_years(ii), numel(alts)), alts, ADEPT_results(:,ii), "magenta.", "DisplayName", "ADEPT")
        else
            scatter3( repelem(ADEPT_years(ii), numel(alts)), alts, ADEPT_results(:,ii), "magenta.", 'HandleVisibility','off')
        end
    end
    hold off
    legend('Location', 'southoutside')
    
end


%% Check Launch Functions

cust_launch_funcs_species = [S_species Sns_species B_species];
MOCAT_SSEM_Launch = zeros(numel(alts), numel(my_sim.results.T));
for i = 1:length(cust_launch_funcs_species) % For Species.   
    cust_species = cust_launch_funcs_species(i);
    cust_name = cust_species.species_properties.sym_name;
    launch_funcs = cust_species.species_properties.lambda_funs;
    figure()
    disp("Species " + cust_name + " " + string(i) + " / " + string(numel(cust_launch_funcs_species)))
    for j = 1:numel(alts) % For each Altitude
        disp("   Shell " + string(j) + " / " + string(numel(alts)))
        cust_species.species_properties.lambda_funs{j}(1, my_sim.results.T);

        % ADEPT FLM Check
        adept_totals = [];
        for m = 2:numel(my_sim.results.T)
            total = [0];
            for n = 1:m-1
                total = total + FLM_steps{n}{j,cust_name};
            end
            adept_totals(m) = total;
        end

        % MOCAT SSEM check.

        alt_time_func = @(t) cust_species.species_properties.lambda_funs{j}(1, t);
        cum_tot = [];
        for l = 1:numel(my_sim.results.T)
            cum_tot(l) = integral(alt_time_func, 1, my_sim.results.T(l));
        end
        if j == 1
            HandleVisibility = "on";
        else
            HandleVisibility = "off";
        end

        hold on
        plot(my_sim.results.T, cum_tot, "DisplayName", "Fit MOCAT-SSEM Launch Function", 'HandleVisibility',HandleVisibility);
        plot(my_sim.results.T(1:end), adept_totals, "DisplayName", "ADEPT Launch Totals Launch Function", 'HandleVisibility',HandleVisibility);
        xlabel('Time (years)');
        ylabel('Cumulative Launch Total');
        title('Cumlative Launch of Species '+strrep(cust_species.species_properties.sym_name + " at altitude " + string(alts(j)) + " km",'p','.'), 'Interpreter', 'none');
        legend('Location', 'southoutside');
        grid on;
        hold off;
    end
end


%MOCAT-SSEM Total
%total = 0;
%for k = 1:numel(FLM_steps)
%    total = total + FLM_steps{k}{j, cust_name};
%end
