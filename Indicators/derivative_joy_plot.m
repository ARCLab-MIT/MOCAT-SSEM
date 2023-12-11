function [figures] = derivative_joy_plot(obj, species, save_folder, min_fit_index, max_fit_index)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    %%% Derivative Joy Plot
    if nargin < 2
        species = obj.scen_properties.species
    end
    
    figures = [];
    N_shell = obj.scen_properties.N_shell;
    for jj = 1:numel(species)
        species_name = species(jj).species_properties.sym_name;
        spec_i = get_species_index_from_name(obj.scen_properties, species_name);
        times = obj.results.T;
        time_mat = repmat(times, [1, N_shell]);
        pop_data = obj.results.species_results_dict(spec_i).value;
        pop_diff = diff(pop_data);
        time_intervals = times(2:end) - times(1:end-1);
        approx_derivatives = pop_diff ./ time_intervals;
        offset = max(approx_derivatives,[],'all');
        %offset = max(offset, 1000);
        %offset = 10 %abs(mean(mean(approx_derivatives)));

        % Create the plot
        if nargin < 4
            min_fit_index = floor(.75*numel(times));
        end
        if nargin < 5
            max_fit_index = numel(times);
        end
        fit_times = times(min_fit_index+1:end);
        
        fit_data = approx_derivatives(min_fit_index:end,:);
        
        % Add a linear fit
        degree = 1; % Degree of the polynomial (1 for linear)]
        fit_time_mat = repmat(fit_times, [1, N_shell]);
        
        coeffs = zeros(N_shell,2);
        fit_values = zeros(N_shell,numel(fit_times));
        for ii = 1:N_shell
            coeffs(ii,:) = polyfit(fit_time_mat(:,ii), fit_data(:,ii), degree);
            fit_values(ii,:) = polyval(coeffs(ii,:) , fit_time_mat(:,ii).');
        end
       
        figures(spec_i) = figure();
        joyPlot(approx_derivatives, obj.results.T(2:end), 0.25*offset, ...
            'reverse', true, 'FaceColor', coeffs(:,1), 'FaceAlpha', 0.6);
        yticklabels(obj.scen_properties.HMid)
        colormap(centered)
        cbar = colorbar();
        cbar.Label.String = "Linear Best Fit Coefficent, T = [150,200]";
        xlabel('Time (years)', 'FontSize', 14);
        ylabel('Altitude Bin Center [km]', 'FontSize', 14);
        cur_ylim = ylim;
        ylim([0, cur_ylim(2)]); %TODO make this variable, % display figures tiled nicely
        title(species_name, 'Interpreter', 'none');
        if exist("save_folder", 'var')
            saveas(gcf, fullfile(save_folder, species_name + " derivative joy plot.fig"))
            saveas(gcf, fullfile(save_folder, species_name + " derivative joy plot.png"))
        end
    end
end