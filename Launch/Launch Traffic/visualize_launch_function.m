function visualize_launch_function(cust_launch_funcs_species, my_sim)
% This function is used to visualize synthetic launch rates created using 
% launch_func_lambda_fun.
% launch_func_lambda_fun is a list of species objects with custom launch
% functions
% my_sim is the simulation.
%% Visualize Launch Rates
%cust_launch_funcs_species = [S_species Su_species Sns_species B_species];
alts = my_sim.scen_properties.HMid;
for i = 1:length(cust_launch_funcs_species) % For Species.   
    cust_species = cust_launch_funcs_species(i);
    cust_name = cust_species.species_properties.sym_name;
    
    launch_data = zeros(numel(my_sim.scen_properties.N_shell), numel(my_sim.scen_properties.scen_times));
    launch_funcs = cust_species.species_properties.lambda_funs;

    disp("Species " + cust_name + " " + string(i) + " / " + string(numel(cust_launch_funcs_species)))
    for j = 1:numel(alts) % For each Altitude
        disp("   Shell " + string(j) + " / " + string(numel(alts)))
        cust_species.species_properties.lambda_funs{j}(1, my_sim.scen_properties.scen_times);
        alt_time_func = @(t) cust_species.species_properties.lambda_funs{j}(1, t);

        if j == 1
            HandleVisibility = "on";
        else
            HandleVisibility = "off";
        end

        launch_data(j,:) =  alt_time_func(my_sim.scen_properties.scen_times);
    end
    figure()
    hold on
    [X, Y] = meshgrid(my_sim.scen_properties.scen_times, alts);
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
end