    clear
    clc
    Cd = 2.2;
    radius = 1.75;
    mass = 1500;
    starting_altitude = 1200;
    a = pi * radius^2;
    deltat = 4;
    pmd = .90;
    TLEs.Alts = [200, 300, 1000];
    points = 365*4;
    sats = 125*4;
  

    scen_properties = MOCATSSEM_Scen_Prop_Cons('start_date', datetime(2023,1,3,0,0,0), ...
                                           'simulation_duration', 200, 'N_steps', ...
                                           200, 'dens_model', @JB2008_dens_func, ...
                                           'integrator', @ode15s, 'N_shell', 40, ...
                                           'h_max', 1400, 'h_min', 200);
    
    candidate = {};
    candidate.adept_id = 9999999999;
    candidate.epoch_start = juliandate(today("datetime"));
    candidate.sma = scen_properties.re + starting_altitude;
    candidate.ecc = 0.0001;
    candidate.inc = 52.0;
    candidate.raan = 360;
    candidate.aop = 120;
    candidate.ma = 45;
    candidate.epoch_end = candidate.epoch_start + 365.25 * deltat;
    candidate.obj_type = 2;
    candidate.disp_option = 0;
    candidate.stationkeeping = 6;
    candidate.area = a;
    candidate.mass = mass;
    candidate.size = radius*2;
    candidate.weight = sats/points;
    candidate = struct2table(candidate);
    candidate_table = repelem(candidate(:,:), points,1);
    candidate_table.weight = candidate_table.weight * (1-pmd);
    candidate_table.epoch_start = candidate_table.epoch_start + 365.25*linspace(0,1,points).';
    candidate_table.epoch_end = candidate_table.epoch_start + 365.25 * deltat;

   
    func = @(x) datetime(x,'convertfrom','juliandate');
    epoch_start_datime = varfun(func,candidate_table,'InputVariables','epoch_start');
    candidate_table.epoch_start_datime_matlab = epoch_start_datime.(1);

    disp("Adding candidate satellites...")
    candidate_cycle_launches = {};
    data_duration = years(days(max(candidate_table.epoch_start) - min(candidate_table.epoch_start))); % Assumes traffic replenished going forward
    repetitions = ceil(scen_properties.simulation_duration/data_duration);
    for repitition = 0:repetitions
        candidate_table_copy = candidate_table;
        start_date = scen_properties.start_date + years(data_duration * repitition);
        end_date = start_date + years(data_duration);
        disp("Regenerating launches for " + string(start_date) + " to " + string(end_date) + ".")
        num_dates = height(candidate_table_copy);
        candidate_table_copy.epoch_start_datime_matlab = candidate_table_copy.epoch_start_datime_matlab + years(data_duration * repitition);
        candidate_table_copy.lifetime = candidate_table_copy.epoch_end - candidate_table_copy.epoch_start;
        candidate_table_copy.epoch_start = juliandate(candidate_table_copy.epoch_start_datime_matlab);
        candidate_table_copy.epoch_end = candidate_table_copy.epoch_start + candidate_table_copy.lifetime;
        candidate_table_copy = removevars(candidate_table_copy, {'epoch_start_datime_matlab', 'lifetime'});
        candidate_cycle_launches{repitition + 1} = candidate_table_copy;
    end

    candidate_cycle_launches= vertcat(candidate_cycle_launches{:});
    candidate_cycle_launches = candidate_cycle_launches(candidate_cycle_launches.epoch_start >= juliandate(scen_properties.start_date) & candidate_cycle_launches.epoch_start <= juliandate(scen_properties.start_date + years(scen_properties.simulation_duration)), :);
    writetable(candidate_cycle_launches, "candidate_cycle_launches.csv")

    %my_sim = applicant_drag_check(scen_properties, mass, radius, starting_altitude, "candidate_cycle_launches.csv");
    
    disposal_alt = 900;
    my_sim = applicant_disposal_alt_check(scen_properties, mass, radius, starting_altitude, disposal_alt, "candidate_cycle_launches.csv");
%%
    % F1
    atmo = " (" + func2str(scen_properties.dens_model) + ")";
    figure1 = my_sim.pop_vs_time_vis();
    title("De-orbit for " + atmo, 'Interpreter', 'none');
    legend("Location", "SouthOutside")
    figure1.Position = [100 100 800 600];

    % F2
    figure2 = my_sim.total_species_evol_vis();
    figure2.Position = [1000 100 800 600];
    view(3); % 3D view

    % F4
    % Calculate the trapz values for each row (population snapshot)
    times = scen_properties.scen_times;
    populationCounts = sum(my_sim.results.X(:,:),2);
    
    % % Create the plot
    figure3 = figure();
    fit_times = times(150:end);
    fit_data = populationCounts(150:end);
    plot(times, populationCounts, 'b', 'Marker', '.');
    xlabel('Time');
    ylabel('Satellites');
    title('Fit', 'Interpreter', 'none');

    % Add a linear fit
    degree = 1; % Degree of the polynomial (1 for linear)
    coeffs = polyfit(fit_times, fit_data, degree);

    % Generate y values for the linear fit
    fit_values = polyval(coeffs, fit_times);

    hold on;
    plot(fit_times, fit_values, 'r', 'LineWidth', 2);

    % Display linear fit equation as text annotation
    fit_eq = sprintf('Linear Fit: y = %.4fx + %.4f', coeffs(1), coeffs(2));
    text(times(end), fit_values(end), fit_eq, 'HorizontalAlignment', 'right');

    legend('N_1500kg', 'Linear Fit', 'Location', 'best', 'Interpreter', 'none');
    hold off;
    figure3.Position = [2000 100 800 600];


    figure4 = figure();
    hold on
    tot_pop = sum(my_sim.results.X(:,:),2);
    times = my_sim.results.T(1:end-1);
    
    % Numerical Derivative
    pop_diff = tot_pop(2:end) - tot_pop(1:end-1);
    time_intervals = my_sim.results.T(2:end) - my_sim.results.T(1:end-1);
    approx_derivatives = pop_diff ./ time_intervals;
    
    % Calculate the percent rate of change
    percent_rate_of_change = (approx_derivatives ./ tot_pop(1:end-1)) * 100;
    plot(times, percent_rate_of_change, 'b', 'Marker', '.');
    yline(5, 'Color', 'red', 'LineWidth', 3, 'DisplayName', "Constraint = 5 %")
    xlabel('Time');
    ylabel('Satellites % change');
    xlim([0 max(my_sim.results.T)])
    ylim([-20 20])
    title('Percent Change', 'Interpreter', 'none');
    hold off
    figure4.Position = [100 800 800 600];



    figure5 = figure();
    hold on
    tot_pop = sum(my_sim.results.X(:,:),2);
    times = my_sim.results.T(1:end-1);
    
    % Numerical Derivative
    pop_diff = tot_pop(2:end) - tot_pop(1:end-1);
    time_intervals = my_sim.results.T(2:end) - my_sim.results.T(1:end-1);
    approx_derivatives = pop_diff ./ time_intervals;
    
    % Calculate the percent rate of change
    percent_rate_of_change = (approx_derivatives ./ tot_pop(1:end-1)) * 100;
    plot(times, approx_derivatives, 'b', 'Marker', '.');
    %yline(5, 'Color', 'red', 'LineWidth', 3, 'DisplayName', "Constraint = 5 %")
    xlabel('Time');
    ylabel('Satellites change');
    title('Derivative', 'Interpreter', 'none');
    xlim([0 max(my_sim.results.T)])
    ylim([-20 20])
    hold off
    figure5.Position = [1000 800 800 600];


        % % F3
    % figure3 = figure();
    % % Calculate the time step
    % timeStep = mean(diff(scen_properties.scen_times));
    % populationCounts = my_sim.results.X(:,:);
    % 
    % % Calculate the numerical derivative of population counts with respect to time
    % populationDerivative = diff(populationCounts, 1, 1) / timeStep;
    % %populationDerivative = 100* populationDerivative./my_sim.results.X(2:end,:);
    % 
    % % Create meshgrid for plotting
    % [altitudeGrid, timeGrid] = meshgrid(altitudes, time(2:end));
    % 
    % % Plot the derivative
    % surf(altitudeGrid, timeGrid, populationDerivative, 'EdgeColor', 'none');
    % xlabel('Altitude');
    % ylabel('Time');
    % zlabel('Population Derivative');
    % title('Numerical Derivative of Population with Respect to Time');
    % colorbar;
    % view(3); % 3D view
    % figure3.Position = [2000 100 800 600];
    % % Adjust axis limits for better visualization
    % xlim([0,800]);
    % ylim([0, 200]);
    % zlim([0, 10]);


    %% Visualize Launch Rates
cust_launch_funcs_species = my_sim.scen_properties.species;
alts = my_sim.scen_properties.HMid;
for i = 1:length(cust_launch_funcs_species) % For Species.   
    cust_species = cust_launch_funcs_species(i);
    cust_name = cust_species.species_properties.sym_name;
    
    launch_data = zeros(numel(my_sim.scen_properties.N_shell), numel(my_sim.results.T));
    launch_funcs = cust_species.species_properties.lambda_funs;

    disp("Species " + cust_name + " " + string(i) + " / " + string(numel(cust_launch_funcs_species)))
    for j = 1:numel(alts) % For each Altitude
        disp("   Shell " + string(j) + " / " + string(numel(alts)))
        cust_species.species_properties.lambda_funs{j}(1, my_sim.results.T);
        alt_time_func = @(t) cust_species.species_properties.lambda_funs{j}(1, t);

        if j == 1
            HandleVisibility = "on";
        else
            HandleVisibility = "off";
        end

        launch_data(j,:) =  alt_time_func(my_sim.results.T);
    end
    figure6 = figure();
    hold on
    [X, Y] = meshgrid(my_sim.results.T, alts);
    surf(X, Y, launch_data,'edgecolor','none');
    xlabel('Time (years)');
    ylabel('Altitude (km)');
    zlabel('Launches Sats/Year');
    title('Launches of Species '+strrep(cust_species.species_properties.sym_name,'p','.'), 'Interpreter', 'none');
    %legend('Location', 'southoutside');
    c = colorbar();
    c.Label.String = "Satellites per Year";
    grid on;
    hold off;
    %view(3)
    figure6.Position = [2000 800 800 600];
end
