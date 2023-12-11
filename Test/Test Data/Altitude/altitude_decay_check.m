
%% Make Scenario Properties
clear
clc
cd(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(fullfile("..", ".."))

VAR = MOCATSSEM_VAR_Cons();

%mass = 1.33; radius = 0.1; starting_altitude = 800; goal_lifetime = 25;
%drag_check(scen_properties, mass, radius, starting_altitude)

%% Test
file = fullfile(".", "TLEs_51470_v4.mat");
mass = 260;
beta = 0.135;
compare_decay(mass, beta, file)


%% Starlink
file = fullfile(".", "TLEs_51470_v4.mat");
mass = 260;
beta = 0.135;
compare_decay(mass, beta, file)

%% Humanity Star
file = fullfile(".", "TLEs_Hstar.mat");
mass = 10.5;
beta = 0.1421;
compare_decay(mass, beta, file)

%% Cubesat
file = fullfile(".", "Cubesat_TLE_JDates_Altitudes_TLE42.mat");
mass = 1.3;
beta = 0.5545;
compare_decay(mass, beta, file)

% opts = detectImportOptions("Cubesat_TLE_JDates_Altitudes_TLE42.csv");
% T = readtable("Cubesat_TLE_JDates_Altitudes_TLE42.csv",opts)
% 
% data = {}
% data.jdate = T{:,"jdate"}
% data.Alts = T{:,"altitude"}
%% 
function my_sim = compare_decay(mass, beta, file)
    
    load(file);
    TLEs.datetimes = datetime(TLEs.jdate,'convertfrom','juliandate');
    TLEs.times = years(TLEs.datetimes - min(TLEs.datetimes));
    start_date =  min(TLEs.datetimes);
    
    starting_altitude = TLEs.Alts(1);
    Cd = 2.2;
    radius = sqrt(Cd*pi*beta*mass)/(Cd * pi);
    a = pi * radius^2;
    
    VAR = MOCATSSEM_VAR_Cons_ALT(TLEs);
    % Chose atmo model
    %VAR.density_filepath = fullfile(".", "Atmosphere Model", "JB2008", "Precomputed", "dens_highvar_2000.mat");
    %VAR.dens_model = @JB2008_dens_func;
    %VAR.time_dep_density = true;
    
    VAR.dens_model = @static_exp_dens_func;
    VAR.time_dep_density = false;
    
    VAR.integrator = @ode15s; % Fixed options: ode1, ode2, ode3, ode4, ode5
    %                                     % Variable options: ode45, ode21,
    %                                     % ode113, oder78, ode89, ode15s,
    %                                     % oder23s, ode23t, ode23tb
    VAR.integrator_type = 'variable'; %or 'fixed';
    
    % Add the time stuff needed
    VAR.start_year = year(start_date);
    VAR.start_month = month(start_date);
    VAR.start_day = day(start_date);
    VAR.start_hour = hour(start_date);
    VAR.start_minute = minute(start_date);
    VAR.start_second = second(start_date);
    VAR.N_steps = 400;
    
    % Get re-entry time from data
    [val,idx]=min(abs(TLEs.Alts-VAR.h_min));
    reentry_time=TLEs.times(idx);
    VAR.simulation_duration = 5*reentry_time; %Years
    VAR.scen_times = linspace(0,VAR.simulation_duration,VAR.N_steps)';
    
    [folder, filename] = fileparts(file);

    % Run Sim
    scen_properties = scen_properties_class(VAR);
    my_sim = drag_check(scen_properties, mass, radius, starting_altitude);

    atmo = " (" + func2str(scen_properties.dens_model) + ")";

    figure = my_sim.pop_vs_time_vis();
    xline(reentry_time, "LineWidth", 5, "LineStyle", ":", "DisplayName", "Actual Re-Entry Time")
    title("De-orbit for " + filename + atmo, 'Interpreter', 'none');
    legend("Location", "SouthOutside")
    figure.Position = [100 100 800 600];
    saveas(figure, "total_pop" + filename + atmo + ".png")
    saveas(figure, "total_pop" + filename + atmo + ".fig")
    

    figure2 = my_sim.total_species_evol_vis();
    ylim([scen_properties.h_min 1.1* starting_altitude])
    xline(reentry_time, "LineWidth", 5, "LineStyle", ":", "DisplayName", "Actual Re-Entry Time")
    legend("Location", "SouthOutside")
    title("De-orbit for " + filename + atmo, 'Interpreter', 'none');
    figure2.Position = [100 100 800 600];
    saveas(figure2, "density_" + filename + atmo + ".png")
    saveas(figure2, "density_" + filename + atmo + ".fig")
    % Set Time one.

end