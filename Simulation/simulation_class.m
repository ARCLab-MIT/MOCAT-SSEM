classdef simulation_class < handle
    properties
        species_list
        scen_properties
        equations
        results
        xdot_eqs % the equations. Used during optimization
        xdot_fun %Used for optimization of run-time by pre-computing non-time-varying equations.
        var_col % Column of state variables for the full system
        drag_term_upper % the drag term corresponding to values entering a shell from above (without density [function]
        drag_term_cur % the drag term correspoding to values leaving a shell to below (without density) [function]
    end % end properties
    methods
        function obj = simulation_class(species_list, scen_properties)
            obj.species_list = species_list;
            obj.scen_properties = scen_properties;
            obj.equations = zeros(scen_properties.N_shell, width(obj.species_list), 'sym');
        end

        function obj = replan_time(obj)
            % This function adjusts aspects of the class that depend on the
            % start_year, start_month, start_day, start_hour, start_minute, 
            % start_second, start_second_fraction, start_date, 
            % simulation_duration, and N_steps parameters.

            % TODO: currently does not adjust launch splines. Add checking
            % for if current spline boundaries include new timespan,
            % recalculate if not.
            obj.scen_properties.scen_times = linspace(0,obj.scen_properties.simulation_duration,obj.scen_properties.N_steps)';

            %% Density stuff
            if isprop(obj.scen_properties, 'density_filepath') && ~isempty(obj.scen_properties.density_filepath)
                [~, name,~] = fileparts(obj.scen_properties.density_filepath);
                dens_var = load(obj.scen_properties.density_filepath, name).(name);
                %interp1(dens_var.alt.', dens_var.dens(:,1), obj.scen_properties.HMid.', 'nearest', 'extrap');
    
                % Construct an array of datetimes based on dens_var times and
                % then find relative time from start_date.
                datetimes = NaT(1, length(dens_var.year));
                for i = 1:length(dens_var.year)
                    datetimes(1, i) = datetime(dens_var.year(i), dens_var.month(i), 1);
                end
                times = years(datetimes - obj.scen_properties.start_date);
        
                % Don't extrapolate under alt range, warn over alt range
                % Being able to nearest-neighbor extrapolate would be nice, but
                % does not appear to be supported by interp2
                if min(dens_var.alt) > min(obj.scen_properties.HMid)
                    error("Minimum model shell median altitude " + ...
                        "(min(scen_properties.HMid)) is below the minimum " + ...
                        "altitude included in the file " + ...
                        obj.scen_properties.density_filepath + ". Please select a " + ...
                        "higher scen_properties.h_min or provide a more " + ...
                        "extensive atmospheric density file")
                end
                if max(dens_var.alt) < max(obj.scen_properties.HMid)
                    warning("Maximum model shell median altitude " + ...
                        "(max(scen_properties.HMid)) exceeds the " + ...
                        "altitude included in the file " + ...
                        obj.scen_properties.density_filepath + ". Density of 0 " + ...
                        "is being used for higher altitudes, but " + ...
                        "may effect result accuracy. Please select a " + ...
                        "lower scen_properties.h_max or provide a more " + ...
                        "extensive atmospheric density file.")
                end
    
                obj.scen_properties.interp_dens_var = interp2(times, dens_var.alt.', dens_var.dens, ...
                        obj.scen_properties.scen_times.', obj.scen_properties.R02.', ...
                        'nearest', 0);
            end % End density interpolation
        end

        function obj = build_model(obj)
            disp("Building model...")
            t = sym('t');
            
            %% Density stuff
            replan_time(obj) % Addresses time-dependent density consideration.

            %%
            full_Cdot_PMD = zeros(obj.scen_properties.N_shell, width(obj.species_list), 'sym');
            %full_drag = sym2cell(zeros(obj.scen_properties.N_shell, width(obj.species_list), 'sym'));
            full_lambda = sym2cell(zeros(obj.scen_properties.N_shell, width(obj.species_list), 'sym'));
            full_coll = zeros(obj.scen_properties.N_shell, width(obj.species_list), 'sym');
                    
            % For optimization we make the function for everything we can
            for i = 1:length(obj.species_list) %Species Index
                var(:,i) = obj.species_list(i).species_properties.sym.';
                num_eqs = size(var(:,i), 1);
                obj.species_list(i).species_properties.eq_idxs = 1 + (num_eqs) * (i-1):(i)*num_eqs;
            end
            var = reshape(var, [obj.scen_properties.N_shell*width(obj.scen_properties.species),1]);
            %TODO: Switch all scenario.var_col to scenario.scen_properties.var_col
            %obj.var_col = var;
            obj.scen_properties.var_col = var;

            % Drag is handled in a bit of a confusing way for computational
            % efficency.  The symbolic terms are computed in this function.
            % If a non-time-varying atmosphere model is used, they are just
            % added to the bigger set of equations and converted to a
            % function at the end.  If a time varying atmosphere model is
            % used, then they are turned into functions and stored as
            % object properties to be multiplied by the time-varying
            % density at run-time.
            drag_term_upper = zeros(obj.scen_properties.N_shell, width(obj.species_list), 'sym');
            drag_term_cur = zeros(obj.scen_properties.N_shell, width(obj.species_list), 'sym');
            
            for i = 1:length(obj.species_list) %Species Index
                disp("Now processing " + obj.species_list(i).species_properties.sym_name)
               
                % Lambda 
                lambda = obj.species_list(i).launch_func(obj.scen_properties.scen_times, ...
                                                         obj.scen_properties.HMid, ...
                                                         obj.species_list(i).species_properties, ...
                                                         obj.scen_properties);
                full_lambda(:,i) = lambda;
              
                %  Post-Mission Disposal
                Cdot_PMD = obj.species_list(i).pmd_func(t, obj.scen_properties.HMid, obj.species_list(i).species_properties, obj.scen_properties);
                % obj.equations(:,i) = obj.equations(:,i) + C_PMD;
                full_Cdot_PMD(:,i) = Cdot_PMD;

                % Drag
                [upper_term, current_term] = drag_func_sym(obj.species_list(i).species_properties, obj.scen_properties);
                try
                    drag_term_upper(:,i) = upper_term;
                    drag_term_cur(:,i) = current_term;
                catch
                    disp("why")
                end

                %if obj.scen_properties.time_dep_density == false
                %    Fdot = obj.species_list(i).drag_func(t, obj.species_list(i).species_properties, obj.scen_properties);
                %    full_drag(:,i) = Fdot;
                %end
                %obj.equations[spec_index] = obj.equations[spec_index] + lambda

            end% end species iteration

            % Collisions
            % full_coll = zeros(obj.scen_properties.N_shell, width(obj.scen_properties.species), 'sym'); %Rows are heights, cols are species.
            for i = 1:length(obj.scen_properties.species_pairs)
                full_coll = full_coll + obj.scen_properties.species_pairs(i).eqs;
            end
            
            obj.equations = full_Cdot_PMD + full_coll;
            
            %Figure out if lambda can be sym equation array, or needs to be
            %cell.
            if IsTypeConsistent(full_lambda, "sym")
                obj.equations = obj.equations + full_lambda;
                obj.scen_properties.sym_lambda = true;
            else
                obj.scen_properties.sym_lambda = false;
                full_lambda = reshape(full_lambda, [obj.scen_properties.N_shell*width(obj.scen_properties.species),1]);
                for i = 1:length(full_lambda)
                    %disp(class(full_lambda{i,1}))
                    if class(full_lambda{i,1}) == "sym"
                        try
                            full_lambda{i,1} = double(full_lambda{i,1});
                        catch
                            disp("Unable to coerce " + string(full_lambda{i,1}) ...
                                  + " to double. This is expected for " + ...
                                  "symbolic launch rates, but results " + ...
                                  "in slow-down.")
                        end %convert to double if possible.
                    end %for sym values
                end %end full lambda loop
                obj.scen_properties.full_lambda = full_lambda;
            end
            % Take any non-function values and make into dummy functions.
            % This allows full cell to be evaluated in population_shell.
            if obj.scen_properties.sym_lambda == false
                for i = 1:length(full_lambda)
                    %disp(class(obj.scen_properties.full_lambda{i,1}))
                    if class(obj.scen_properties.full_lambda{i,1}) ~= "function_handle"
                        try
                            val =  obj.scen_properties.full_lambda{i,1};
                            obj.scen_properties.full_lambda{i,1} = @(x, t) val;
                        catch
                            disp(i)
                        end
                    end
                end
            end


%             % Leaves out drag if time dependent, then computed during
%             % integration.
%             if IsTypeConsistent(full_drag, "sym")
%                 obj.equations = obj.equations + full_drag;
%                 obj.scen_properties.sym_drag = true;
%             else
%                 obj.scen_properties.sym_drag = false;
%                 full_drag = reshape(full_drag, [obj.scen_properties.N_shell*width(obj.scen_properties.species),1]);
%                 obj.scen_properties.full_drag = full_drag;
%             end

            if obj.scen_properties.time_dep_density == false
                rho = obj.scen_properties.dens_model(0, obj.scen_properties.R02, obj.scen_properties);
                rho_mat = repelem(rho, width(obj.scen_properties.species));
                rho_mat = reshape(rho_mat.', [width(obj.scen_properties.species), obj.scen_properties.N_shell + 1]).';
                full_drag = drag_term_upper .* rho_mat(2:end,:) + drag_term_cur .* rho_mat(1:end-1, :);
                obj.scen_properties.full_drag = full_drag;
                obj.equations = obj.equations + full_drag;
                obj.scen_properties.sym_drag = true;
                
            elseif obj.scen_properties.time_dep_density == true
                 drag_term_upper = reshape(drag_term_upper, [obj.scen_properties.N_shell*width(obj.scen_properties.species),1]);
                 drag_term_cur = reshape(drag_term_cur, [obj.scen_properties.N_shell*width(obj.scen_properties.species),1]);
                 obj.drag_term_upper = matlabFunction(drag_term_upper,'Vars',{obj.scen_properties.var_col });
                 obj.drag_term_cur =  matlabFunction(drag_term_cur,'Vars',{obj.scen_properties.var_col });
            else
                error("Please set property scen_properties.time_dep_density " + ...
                    "to true if using JB2008 or another time-dependent density model," + ...
                    "or false if using a static atmospheric density model.")
            end
    
            % Turn obj.equations into a matlabFunction for speed
            obj.xdot_eqs =  reshape(obj.equations, [obj.scen_properties.N_shell*width(obj.scen_properties.species),1]);


            % Make Integrated Indicator Variables if passed.
            %TODO Check
            if isprop(obj.scen_properties, 'integrated_indicator_var_list')
                % Parse into equations
                integrated_indicator_var_list = obj.scen_properties.integrated_indicator_var_list;
                for ind_i = 1:length(integrated_indicator_var_list)
                    if isempty(integrated_indicator_var_list(ind_i).eqs)  %This checking is also done in make_indicator_eqs
                        integrated_indicator_var_list(ind_i) = make_indicator_eqs(obj, integrated_indicator_var_list(ind_i));
                    end
                end
                % Assign number indicators as property for number of
                % variables & indexes
                obj.scen_properties.num_integrated_indicator_vars = 0;
                end_indicator_idxs = size(obj.xdot_eqs,1);
                for ind_i = 1:length(integrated_indicator_var_list)
                    % Keep track of how many indicator vars we are adding
                    num_add_indicator_vars = size(integrated_indicator_var_list(ind_i).eqs, 1);
                    obj.scen_properties.num_integrated_indicator_vars = obj.scen_properties.num_integrated_indicator_vars + num_add_indicator_vars;
                    
                    % Keep Track of their indexes in obj.xdot_eqs
                    start_indicator_idxs = end_indicator_idxs + 1;
                    end_indicator_idxs = start_indicator_idxs + num_add_indicator_vars - 1;
                    integrated_indicator_var_list(ind_i).indicator_idxs = (start_indicator_idxs:end_indicator_idxs);
                    
                    % Glue variables onto end of obj.xdot_eqs
                    obj.xdot_eqs = [obj.xdot_eqs; integrated_indicator_var_list(ind_i).eqs];
                end

                % Pad the equations for real-time evaluations, if needed
                % if obj.scen_properties.time_dep_density == true
                %      indicator_pad = zeros(obj.scen_properties.num_integrated_indicator_vars, 1, 'sym');
                %      drag_term_upper = reshape(drag_term_upper, [obj.scen_properties.N_shell*width(obj.scen_properties.species),1]);
                %      drag_term_upper = [drag_term_upper;indicator_pad];
                %      drag_term_cur = reshape(drag_term_cur, [obj.scen_properties.N_shell*width(obj.scen_properties.species),1]);
                %      drag_term_cur = [drag_term_cur;indicator_pad];
                %      obj.drag_term_upper = matlabFunction(drag_term_upper,'Vars',{obj.scen_properties.var_col });
                %      obj.drag_term_cur =  matlabFunction(drag_term_cur,'Vars',{obj.scen_properties.var_col });
                % end
                if obj.scen_properties.sym_lambda == false
                    indicator_pad = cell(obj.scen_properties.num_integrated_indicator_vars, 1);
                    indicator_pad(:) = {@(x, t) 0};
                    obj.scen_properties.full_lambda = [obj.scen_properties.full_lambda;indicator_pad];
                end

            end

            % Make Non-Integrated Indicator Variables if passed.
            % No need for all the dimension adjustments for indicators
            % integrated as part of xdot.
            %TODO Check
            if isprop(obj.scen_properties, 'indicator_var_list')
                % Parse into equations
                indicator_var_list = obj.scen_properties.indicator_var_list;
                for ind_i = 1:length(indicator_var_list)
                    if isempty(indicator_var_list(ind_i).eqs)  %This checking is also done in make_indicator_eqs
                        indicator_var_list(ind_i) = make_indicator_eqs(obj, indicator_var_list(ind_i));
                    end
                end
            end

            obj.xdot_fun = matlabFunction(obj.xdot_eqs,'Vars',{obj.scen_properties.var_col });

            
            disp("Done building model.")
        end %end build model

        function results = run_model(obj, x0, varargin) 
            %progressBar (bool) whether to display progress bar
            % x0 is state as either array or vector, with or without
            % indicator vars/
            
            % TODO: Fix input parser so that fields manually set to
            % scen_properties default to set property not generic default
            % below.
            % start_date (datetime) simulate start time
            % simulation_duration (double) simulation duration in years
            % N_step (double) number of steps to show results for in
            % simulation range.
            % disp_times (bool) show the time steps to console during
            % integration.
            p = inputParser;
            addRequired(p, 'obj', @(x)class(x) == "simulation_class")
            state_elms = obj.scen_properties.N_shell * size(obj.scen_properties.species, 2);
            addRequired(p, 'x0', @(x) numel(x0) == state_elms | numel(x0) == state_elms + obj.scen_properties.num_integrated_indicator_vars)
            addOptional(p, 'progressBar', true, @islogical)
            addOptional(p, 'disp_times', true, @islogical)
            addOptional(p, 'start_date', obj.scen_properties.start_date, @isdatetime)
            addOptional(p, 'simulation_duration', obj.scen_properties.simulation_duration, @(x)class(x) == "double")
            addOptional(p, 'N_step', obj.scen_properties.N_step, @(x) isinteger(uint16(x)) & x>0)
            parse(p,obj, x0, varargin{:});
          
            adjustedTime = false;

            if p.Results.start_date ~= obj.scen_properties.start_date
                obj.scen_properties.start_year = year(p.Results.start_date);
                obj.scen_properties.start_month = month(p.Results.start_date);
                obj.scen_properties.start_day = day(p.Results.start_date);
                obj.scen_properties.start_hour = hour(p.Results.start_date);
                obj.scen_properties.start_minute = minute(p.Results.start_date);
                obj.scen_properties.start_second = floor(second(p.Results.start_date));
                obj.scen_properties.start_second_fraction = mod(second(p.Results.start_date), 1);
                obj.scen_properties.start_date = p.Results.start_date;
                adjustedTime = true;
            end

            if p.Results.simulation_duration ~= obj.scen_properties.simulation_duration
                obj.scen_properties.simulation_duration = p.Results.simulation_duration;
                adjustedTime = true;
            end

           if p.Results.N_step ~= obj.scen_properties.N_step
                obj.scen_properties.N_step = p.Results.N_step;
                adjustedTime = true;
           end
            
           if adjustedTime
               disp("Replanning time...")
               obj.replan_time();
           end

           if obj.scen_properties.disp_times
               disp("Running model...")
           end
    
            % Populate matrix with all the variables 
            % Resulting matrix is N_shell x size(species_list)
            for i = 1:length(obj.species_list) %Species Index
                var(:,i) = obj.species_list(i).species_properties.sym.';
            end

            %var = reshape(var, [obj.scen_properties.N_shell*width(obj.scen_properties.species),1]);
            %int_equations =  reshape(obj.equations, [obj.scen_properties.N_shell*width(obj.scen_properties.species),1]);
            %fun = matlabFunction(int_equations,'Vars',{var});
            %func=@(t,x) fun(x);
            %func = population_shell(time,x,scen_properties, int_equations, var)

            % Progress bar
            if p.Results.progressBar
                options_ode = odeset('reltol', 1.e-4,'abstol', 1.e-4, 'OutputFcn',@odeprog,'Events',@odeabort);
            else
                options_ode = odeset('reltol', 1.e-4,'abstol', 1.e-4);
            end
            
            % Display time steps in console or not.
            obj.scen_properties.disp_times = p.Results.disp_times;


            xdot_func = @(x,t)population_shell(x,t,obj);

            % Ravel it if it has the right number
            if numel(x0) == obj.scen_properties.N_shell * numel(obj.scen_properties.species)
                x0 = reshape(x0, [], 1); % convert table to vector
            end

            % Pad with the indicator vars if needed.
            if isprop(obj.scen_properties, 'integrated_indicator_var_list') && numel(x0) == obj.scen_properties.N_shell * numel(obj.scen_properties.species)
                % Pad initial conditions with the indicator variables
                x0 = [x0; zeros(obj.scen_properties.num_integrated_indicator_vars, 1)];
            end
            if obj.scen_properties.disp_times
               tic
            end
            if ~isprop(obj.scen_properties, 'integrator_type')
                 error("Simulation scen_properties.integrator_type must be " + ...
                "set to either 'variable' for variable timestep " + ...
                "integrators (e.g. ode45) or 'fixed' for fixed time step" + ...
                "integrators (e.g. ode5).")
            elseif strcmp(obj.scen_properties.integrator_type,'variable')
                [T,X] = obj.scen_properties.integrator(xdot_func, obj.scen_properties.scen_times, x0, options_ode);
            elseif strcmp(obj.scen_properties.integrator_type,'fixed')
                error("Fixed timestep integrators are not successfully implemented yet.")
                %[X] = obj.scen_properties.integrator(xdot_func, obj.scen_properties.scen_times, x0);
                %[T] = obj.scen_properties.scen_times;
            else
                error("Simulation scen_properties.integrator_type must be " + ...
                    "set to either 'variable' for variable timestep " + ...
                    "integrators (e.g. ode45) or 'fixed' for fixed time step" + ...
                    "integrators (e.g. ode5).")
            end
            if obj.scen_properties.disp_times
               toc
            end
            %TODO: Fix the outputs to be nice.
            obj.results.T = T;
            
            % Split integrated indicator variables from states
            if isprop(obj.scen_properties, 'integrated_indicator_var_list')
                obj.results.integrated_indicators = struct;
                for ind_i = 1:length(obj.scen_properties.integrated_indicator_var_list)
                    integrated_indicator_var = obj.scen_properties.integrated_indicator_var_list(ind_i);
                    obj.results.integrated_indicators.(integrated_indicator_var.name) = X(:, integrated_indicator_var.indicator_idxs);
                end
                obj.results.integrated_indicators_full = X(:,size(X,2)-obj.scen_properties.num_integrated_indicator_vars + 1:end);
                obj.results.X = X(:, 1:size(X,2)-obj.scen_properties.num_integrated_indicator_vars);
            else
                obj.results.X = X;
            end

            % Evaluate non-indicator variables using states
            if isprop(obj.scen_properties, 'indicator_var_list')
                obj.results.indicators = struct;
                for ind_i = 1:length(obj.scen_properties.indicator_var_list)
                    indicator_var = obj.scen_properties.indicator_var_list(ind_i);
                    % Evaluate Indicator Var from results
                    indicator_fun = matlabFunction(indicator_var.eqs,'Vars',{obj.scen_properties.var_col });
                    num_times = size(obj.results.X,1);
                    evaluated_indicator_var = splitapply(indicator_fun,obj.results.X.',1:num_times);
                    obj.results.indicators.(indicator_var.name) = evaluated_indicator_var;
                end
            end
            
            species_array = reshape(obj.results.X,[height(obj.results.X), ...
                    width(obj.results.X)/length(obj.species_list), length(obj.species_list)]);
            species_results_dict = struct;
            for i = 1:length(obj.species_list)
                name = obj.species_list(i).species_properties.sym_name;
                species_results_dict(i).key = name;
                species_results_dict(i).value = species_array(:,:,i);
            end
            obj.results.species_array = species_array;
            obj.results.species_results_dict = species_results_dict;
        end

        function visuals = total_species_evol_vis(obj)
            %Produces a Figure with subplots for each species population 
            % evolution over time as given by the results function
            if isprop(obj,"results")
                species_array = obj.results.species_array;
                alts = obj.scen_properties.HMid;
                disp("Producing visuals for the evolution of each species.")
                if length(obj.species_list)<7
                    figure;
                    pos = get(gcf, 'Position');
                    pos(3) = 600; % change the width
                    pos(4) = 150*length(obj.species_list); % change the height
                    set(gcf, 'Position', pos);
                    for i = 1:length(obj.species_list)
                        %if number of species is less than 7 use subplots
                        %to plot all species evolutions in one figure
                        subplot(length(obj.species_list),1,i)
                        surf(obj.results.T,alts,transpose(species_array(:,:,i)), 'edgecolor','none', "DisplayName", obj.species_list(i).species_properties.sym_name);
                        c = colorbar;
                        c.Label.String = 'Number of Objects';
                        zlabel('Number of Objects')
                        xlabel('Time (years)')
                        ylabel('Altitude (km)')
                        xlim([0 max(obj.results.T)])
                        ylim([min(alts) max(alts)])
                        title('Population of Species '+strrep(obj.species_list(i).species_properties.sym_name,'p','.'), 'Interpreter', 'none')
                        grid on
                        view (0, 90)
                    end
                    visuals = gcf;
                else
                    for i = 1:length(obj.species_list)
                        figure()
                        %if more than 7 species then create visual for 
                        %each species in separate figure
                        surf(obj.results.T,alts,transpose(species_array(:,:,i)),'edgecolor','none', "DisplayName", obj.species_list(i).species_properties.sym_name);
                        c = colorbar;
                        c.Label.String = 'Number of Objects';
                        zlabel('Number of Objects')
                        xlabel('Time (years)')
                        ylabel('Altitude (km)')
                        xlim([0 max(obj.results.T)])
                        ylim([min(alts) max(alts)])
                        title('Population of Species '+strrep(obj.species_list(i).species_properties.sym_name,'p','.'), 'Interpreter', 'none')
                        grid on
                        view (0, 90)
                    end
                end
            else
              error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end 
        end

        function visuals2 = density_species_evol_vis(obj)
            %Produces a Figure with subplots for each species population 
            % DENSITY evolution over time as given by the results function
            if isprop(obj,"results")
                species_array = obj.results.species_array;
                alts = obj.scen_properties.HMid;
                disp("Producing visuals for the evolution of each species' density.")
                if length(obj.species_list)<7
                    figure;
                    pos = get(gcf, 'Position');
                    pos(3) = 600; % change the width
                    pos(4) = 150*length(obj.species_list); % change the height
                    set(gcf, 'Position', pos);
                    for i = 1:length(obj.species_list)
                        %if number of species is less than 7 use subplots
                        %to plot all species evolutions in one figure
                        subplot(length(obj.species_list),1,i)
                        dens_spec = transpose(species_array(:,:,i))./repmat(reshape(obj.scen_properties.V,[],1),1,size(transpose(species_array(:,:,i)),2));
                        surf(obj.results.T,alts,dens_spec,'edgecolor','none')
                        c = colorbar;
                        c.Label.String = 'Density (obj/km^3)';
                        zlabel('Density (obj/km^3)')
                        xlabel('Time (years)')
                        ylabel('Altitude (km)')
                        xlim([0 max(obj.results.T)])
                        ylim([min(alts) max(alts)])
                        title('Population of Species '+strrep(obj.species_list(i).species_properties.sym_name,'p','.'), 'Interpreter', 'none')
                        grid on
                        view (0, 90)
                    end
                    visuals = gcf;
                else
                    for i = 1:length(obj.species_list)
                        figure()
                        %if more than 7 species then create visual for 
                        %each species in separate figure
                        dens_spec = transpose(species_array(:,:,i))./repmat(reshape(obj.scen_properties.V,[],1),1,size(transpose(species_array(:,:,i)),2));
                        surf(obj.results.T,alts,dens_spec,'edgecolor','none')
                        c = colorbar;
                        c.Label.String = 'Density (obj/km^3)';
                        zlabel('Density (obj/km^3)')
                        xlabel('Time (years)')
                        ylabel('Altitude (km)')
                        xlim([0 max(obj.results.T)])
                        ylim([min(alts) max(alts)])
                        title('Population of Species '+strrep(obj.species_list(i).species_properties.sym_name,'p','.'), 'Interpreter', 'none')
                        grid on
                        view (0, 90)
                    end
                end
            else
              error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end 
        end

        function total_pop_time = pop_vs_time_vis(obj)
            rgb_c = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],...
                [0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840]};
            %Produces a Figure for the evolution of the total population of 
            % each species over time across all shells 
            if isprop(obj,"results")
                species_array = obj.results.species_array;
                disp("Producing visual of the change in populations across all altitudes vs time.")
                total_pop_time = figure();
                movegui('center')
                pos = get(gcf, 'Position');
                pos(3) = 600; % change the width
                pos(4) = 800; % change the height
                set(gcf, 'Position', pos);
                for i = 1:length(obj.species_list)
                    if i>length(rgb_c)
                        colorset = [rand rand rand];
                    else
                        colorset = rgb_c{i};
                    end
                    plot(obj.results.T,sum((species_array(:,:,i)),2), LineWidth=2.0, ...
                        DisplayName=strrep(obj.species_list(i).species_properties.sym_name, 'p', '.'),Color=colorset)
                    hold on
                    xlabel('Time (years)')
                    ylabel('Number of Objects')
                    xlim([0 max(obj.results.T)])
                    title('Evolution of Species Population')
                    legend('Location','bestoutside', 'Interpreter', 'none')
                    grid on
                    %set(gca, 'Yscale', 'log')
                end
            else
              error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end 
        end

        function visuals = total_deb_pop_time_deriv(obj, varargin)
            % Graphs the rate of change for the total inactive population 
            % added across all altitudes at each time.
            % percentage = (bool, defaults to false) that expresses rate of
            % change in percentage terms. 
            
            p = inputParser;
            addRequired(p, 'obj', @(x)class(x) == "simulation_class")
            addOptional(p, 'percentage', false, @islogical)
            addOptional(p, 'constraint', 0)
            parse(p,obj, varargin{:});
            percentage = p.Results.percentage;
            percentage_string = "";
            if percentage
                percentage_string = " (%)";
            end
            constraint = p.Results.constraint;
            rgb_c = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],...
                [0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840]};
            %Produces a Figure for the evolution of the total population of 
            % each species over time across all shells 
            if isprop(obj,"results")
                species_array = obj.results.species_array;
                disp("Producing visual of the change in inactive populations across all altitudes vs time.")
                figure()
                movegui('center')
                pos = get(gcf, 'Position');
                pos(3) = 600; % change the width
                pos(4) = 800; % change the height
                set(gcf, 'Position', pos);
                inactive_i = [];
                for i = 1:length(obj.species_list)                    
                    if obj.species_list(i).species_properties.active == false
                        inactive_i = [inactive_i i];
                    end
                end
                for ii = 1:length(inactive_i)
                    i = inactive_i(ii);
                    if i>length(rgb_c)
                        colorset = [rand rand rand];
                    else
                        colorset = rgb_c{i};
                    end
                    disp(obj.species_list(i).species_properties.sym_name)
                    tot_pop = sum((species_array(:,:,i)),2);
                    if percentage
                        diff_tot_pop = 100 * diff(tot_pop)./diff(obj.results.T)./tot_pop(1:end-1,:);
                    else 
                        diff_tot_pop = diff(tot_pop)./diff(obj.results.T);
                    end
                    mid_times = (obj.results.T(1:end-1,:) + obj.results.T(2:end,:)) / 2;
                    plot(mid_times,diff_tot_pop, LineWidth=2.0, ...
                        DisplayName=strrep(obj.species_list(i).species_properties.sym_name, 'p', '.'),Color=colorset)
                    hold on
                    xlabel('Time (years)')
                    ylabel('Change in Number of Objects' + percentage_string)
                    %xlim([0 max(obj.results.T)])
                    %ylim([-20 20])
                    title('Evolution of Debris Population Rate of Change (Increases Only)' + percentage_string)
                    legend('Location','bestoutside', 'Interpreter', 'none')
                    grid on
                    set(gca, 'Yscale', 'log')
                end
                if constraint
                    yline(constraint, 'Color', 'red', 'LineWidth', 3, 'DisplayName', "Constraint = " + constraint + percentage_string)
                end
            else
              error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end 
        end
        function visuals = total_deb_pop_time_deriv2(obj, varargin)
            % Graphs the rate of change for the total inactive population 
            % added across all altitudes at each time.
            % percentage = (bool, defaults to false) that expresses rate of
            % change in percentage terms. 
            
            p = inputParser;
            addRequired(p, 'obj', @(x)class(x) == "simulation_class")
            addOptional(p, 'percentage', false, @islogical)
            addOptional(p, 'constraint', 0)
            parse(p,obj, varargin{:});
            percentage = p.Results.percentage;
            percentage_string = "";
            if percentage
                percentage_string = " (%)";
            end
            constraint = p.Results.constraint;
            rgb_c = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],...
                [0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840]};
            %Produces a Figure for the evolution of the total population of 
            % each species over time across all shells 
            if isprop(obj,"results")
                species_array = obj.results.species_array;
                disp("Producing visual of the change in inactive populations across all altitudes vs time.")
                figure()
                movegui('center')
                pos = get(gcf, 'Position');
                pos(3) = 600; % change the width
                pos(4) = 800; % change the height
                set(gcf, 'Position', pos);
                inactive_i = [];
                for i = 1:length(obj.species_list)                    
                    if obj.species_list(i).species_properties.active == false
                        inactive_i = [inactive_i i];
                    end
                end
                for ii = 1:length(inactive_i)
                    i = inactive_i(ii);
                    if i>length(rgb_c)
                        colorset = [rand rand rand];
                    else
                        colorset = rgb_c{i};
                    end
                    disp(obj.species_list(i).species_properties.sym_name)
                    
                    tot_pop = sum((species_array(:,:,i)),2);
                    times = obj.results.T(1:end-1);
                    
                    % Numerical Derivative
                    pop_diff = tot_pop(2:end) - tot_pop(1:end-1);
                    time_intervals = obj.results.T(2:end) - obj.results.T(1:end-1);
                    approx_derivatives = pop_diff ./ time_intervals;
                    
                    % Calculate the percent rate of change
                    percent_rate_of_change = (approx_derivatives ./ tot_pop(1:end-1)) * 100;
                    
                    % Recover Data from PROC
                    fraction_rate_of_changes = percent_rate_of_change / 100;
                    integrated_values = cumtrapz(fraction_rate_of_changes.*tot_pop(1:end-1));

                    hold on
                    if percentage
                        plot(times,percent_rate_of_change, LineWidth=2.0, ...
                        DisplayName=strrep(obj.species_list(i).species_properties.sym_name, 'p', '.'),Color=colorset)
                    else
                        plot(times,approx_derivatives, LineWidth=2.0, ...
                            DisplayName=strrep(obj.species_list(i).species_properties.sym_name, 'p', '.'),Color=colorset)
                    end
                  
                    % Debug
                    %plot(times,tot_pop(1:end-1), LineWidth=1.0,...
                    %     DisplayName=obj.species_list(i).species_properties.sym_name + " Total Population",Color="k")
                    %plot(times,cumtrapz(approx_derivatives), LineWidth=2.0, LineStyle="--", ...
                    %    DisplayName=strrep(obj.species_list(i).species_properties.sym_name + " Integrated Derivative", 'p', 'o'),Color="b")
                    %plot(times,integrated_values, LineWidth=2.0, LineStyle="--", ...
                    %    DisplayName=strrep(obj.species_list(i).species_properties.sym_name + " Integrated Percentage", 'p', 'x'),Color="r")

                    xlabel('Time (years)')
                    ylabel('Change in Number of Objects' + percentage_string)
                    %xlim([0 max(obj.results.T)])
                    %ylim([-20 20])
                    title('Evolution of Debris Population Rate of Change (Increases Only)' + percentage_string)
                    legend('Location','bestoutside', 'Interpreter', 'none')
                    grid on
                    set(gca, 'Yscale', 'log')
                    hold off
                end
                if constraint
                    yline(constraint, 'Color', 'red', 'LineWidth', 3, 'DisplayName', "Constraint = " + constraint + percentage_string)
                end
            else
              error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end 
        end
        function visuals = total_deb_species_deriv_evolv_vis(obj, varargin)
            %Produces a Figure with subplots for numerical derivative of each 
            % inactive species population evolution over time
            % constraint is the value used in the color scale.
            % percentage is bool (default false) for whether to use percentage or absolute.
            % color_scale_crop_percentiles is pair (default [0, 100]). Color scale will be 
            % set based on the inclusive range from the lower value to the
            % higher value.
            p = inputParser;
            addRequired(p, 'obj', @(x)class(x) == "simulation_class")
            addOptional(p, 'percentage', false, @islogical)
            addOptional(p, 'constraint', 0)
            addOptional(p, 'color_scale_crop_percentiles', [0, 100])
            parse(p,obj, varargin{:});
            percentage = p.Results.percentage;
            constraint  = p.Results.constraint;
            color_scale_crop_percentiles = p.Results.color_scale_crop_percentiles;

            if isprop(obj,"results")
                species_array = obj.results.species_array;
                alts = obj.scen_properties.HMid;
                disp("Producing visuals for the rate change of each debris species.")
                inactive_i = [];
                for i = 1:length(obj.species_list)                    
                    if obj.species_list(i).species_properties.active == false
                        inactive_i = [inactive_i i];
                    end
                end
                mid_times = (obj.results.T(1:end-1,:) + obj.results.T(2:end,:)) / 2;
                if length(inactive_i)<7
                    figure;
                    pos = get(gcf, 'Position');
                    pos(3) = 600; % change the width
                    pos(4) = 150*length(inactive_i); % change the height
                    set(gcf, 'Position', pos);
                    for ii = 1:length(inactive_i)
                        figure;
                        i = inactive_i(ii);
                        if percentage == false
                            diff_pop_per_shell = diff(species_array(:,:,i))./diff(obj.results.T);
                        elseif percentage == true
                            diff_pop_per_shell = 100 * diff(species_array(:,:,i))./diff(obj.results.T)./species_array(1:end-1,:,i);
                        end
                        %if more than 7 species then create visual for 
                        %each species in separate figure
                        surf(mid_times,alts,diff_pop_per_shell.','edgecolor','none')
                        z = diff_pop_per_shell;
                        cmap = flipud(lbmap(256,'BrownBlue'));

                        zlims = prctile( diff_pop_per_shell(:) , color_scale_crop_percentiles);
                        zmin = zlims(1);
                        zmax = zlims(2);

                        zdiff_max = abs(constraint - zmax);
                        zdiff_min = abs(constraint - zmin);
                        clim_bounds = [constraint-max(zdiff_min, zdiff_max), constraint + max(zdiff_min, zdiff_max)];
                        colormap(cmap)
                        c = colorbar();
                        clim(clim_bounds)
                        if percentage == false
                            zlabel('Change in objects per year')
                            c.Label.String = 'Change in objects per year';
                        elseif percentage == true
                            zlabel('Change in objects per year (%)')
                            c.Label.String = 'Change in objects per year (%)';
                        end
                        xlabel('Time (years)')
                        ylabel('Altitude (km)')
                        xlim([0 max(obj.results.T)])
                        ylim([min(alts) max(alts)])
                        title('Species '+strrep(obj.species_list(i).species_properties.sym_name,'p','.'), 'Interpreter', 'none')
                        grid on
                        view (0, 90)
                    end
                    visuals = gcf;
                else
                    for ii = 1:length(inactive_i)
                        i = inactive_i(ii);
                        figure()
                        if percentage == false
                            diff_pop_per_shell = diff(species_array(:,:,i))./diff(obj.results.T);
                        elseif percentage == true
                            diff_pop_per_shell = 100 * diff(species_array(:,:,i))./diff(obj.results.T)./species_array(1:end-1,:,i);
                        end
                        %if more than 7 species then create visual for 
                        %each species in separate figure
                        surf(mid_times,alts,diff_pop_per_shell.','edgecolor','none')
                        z = diff_pop_per_shell;
                        cmap = flipud(lbmap(256,'BrownBlue'));

                        zlims = prctile( diff_pop_per_shell(:) , color_scale_crop_percentiles);
                        zmin = zlims(1);
                        zmax = zlims(2);

                        zdiff_max = abs(constraint - zmax);
                        zdiff_min = abs(constraint - zmin);
                        clim_bounds = [constraint-max(zdiff_min, zdiff_max), constraint + max(zdiff_min, zdiff_max)];
                        colormap(cmap)
                        c = colorbar();
                        clim(clim_bounds)
                        if percentage == false
                            zlabel('Change in objects per year')
                            c.Label.String = 'Change in objects per year';
                        elseif percentage == true
                            zlabel('Change in objects per year (%)')
                            c.Label.String = 'Change in objects per year (%)';
                        end
                        xlabel('Time (years)')
                        ylabel('Altitude (km)')
                        xlim([0 max(obj.results.T)])
                        ylim([min(alts) max(alts)])
                        title('Species '+strrep(obj.species_list(i).species_properties.sym_name,'p','.'), 'Interpreter', 'none')
                        grid on
                        view (0, 90)
                    end
                end
            else
              error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end 
        end

        function total_vs_time = tot_vs_time_vis(obj)
            %Produces a Figure for the evolution of the total population of 
            % all objects over time across all shells 
            if isprop(obj,"results")
                species_array = obj.results.species_array;
                disp("Producing visual of the change in the total population across all altitudes vs time.")
                figure()
                plot(obj.results.T,sum(sum((species_array(:,:,:)),2),3), LineWidth=2.0)
                hold on
                xlabel('Time (years)')
                ylabel('Number of Objects')
                xlim([0 max(obj.results.T)])
                title('Evolution of Total Population')
                grid on
            else
              error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end 
        end

        function species_time = species_time_vis(obj)
            %Produces a Figure per species for the evolution of the total population of 
            % each species over time across all shells 
            rgb_c = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],...
                [0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840]};
            if isprop(obj,"results")
                figure;
                pos = get(gcf, 'Position');
                pos(3) = 1800; % change the width
                pos(4) = 700; % change the height
                set(gcf, 'Position', pos);
                if length(obj.species_list)>5
                    tile_x = round(length(obj.species_list)/4);
                    tile_y = 4;
                    tl = tiledlayout(tile_x,tile_y);
                else
                    tl = tiledlayout(1,length(obj.species_list));
                end
                species_array = obj.results.species_array;
                disp("Producing visuals of the change in each species across all altitudes vs time.")
                for i = 1:length(obj.species_list)
                    if i>length(rgb_c)
                        colorset = [rand rand rand];
                    else
                        colorset = rgb_c{i};
                    end
                    nexttile
                    plot(obj.results.T,sum((species_array(:,:,i)),2), LineWidth=2.0, ...
                        DisplayName=strrep(obj.species_list(i).species_properties.sym_name, 'p', '.'),Color=colorset)
                    xlabel('Time (years)')
                    ylabel('Number of Objects')
                    xlim([0 max(obj.results.T)])
                    title('Evolution of '+strrep(obj.species_list(i).species_properties.sym_name, 'p', '.'),'Interpreter', 'none')
                    grid on
                end
            else
              error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end 
        end
        
        function total_density = density_vs_time(obj)
            % Displays total spatial density in each shell over time 
            if isprop(obj,"results")
                alts = obj.scen_properties.HMid;
                disp("Producing visuals for the evolution of total population density.")
                figure()
                surf(obj.results.T,alts,transpose(sum((obj.results.species_array(:,:,:)),3)./obj.scen_properties.V),'edgecolor','none')
                c = colorbar;
                c.Label.String = 'Density (obj/km^3)';
                zlabel('Density (obj/km^3)')
                xlabel('Time (years)')
                ylabel('Altitude (km)')
                xlim([0 max(obj.results.T)])
                ylim([min(alts) max(alts)])
                title('Total Population Density', 'Interpreter', 'none')
                grid on
                view (0, 90)
            else
                error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end
        end

        function num_time_prod = num_time_prod(obj, x0)
            % Computes and displays Number_Time_Product
            if isprop(obj,"results")
                for spec_ind = 1:length(obj.species_list)
                    spec_launch_save{spec_ind} = obj.species_list(spec_ind).launch_func;
                    obj.species_list(spec_ind).launch_func = @launch_func_null; %changing launch function to 0, note original simulation is unaffected.
                end
                res_array_N = zeros(length(obj.results.T),0);
                res_array_U = zeros(length(obj.results.T),0);
                for i = 1:length(obj.results.species_results_dict)
                        name = obj.results.species_results_dict(i).key;
                        if contains(name,'N')
                            res_array_N = cat(2,res_array_N,sum(obj.results.species_results_dict(i).value,2));
                        elseif contains(name,'U')
                            res_array_U = cat(2,res_array_U,sum(obj.results.species_results_dict(i).value,2));
                        end
                end
    
                my_sim2 = simulation_class(obj.species_list, obj.scen_properties);
                disp('Building simulation for no future launch case.')
                my_sim2.build_model();
                my_sim2.run_model(x0, true);
    %             my_sim2.total_species_evol_vis() %shows the evolution for no future launch case
                for spec_ind = 1:length(obj.species_list)
                    obj.species_list(spec_ind).launch_func = spec_launch_save{spec_ind}; %changing launch function back to original, note original simulation is unaffected.
                end
                res_array_N2 = zeros(length(obj.results.T),0);
                res_array_U2 = zeros(length(obj.results.T),0);
                for i = 1:length(my_sim2.results.species_results_dict)
                        name = my_sim2.results.species_results_dict(i).key;
                        if contains(name,'N')
                            res_array_N2 = cat(2,res_array_N2,sum(my_sim2.results.species_results_dict(i).value,2));
                        elseif contains(name,'U')
                            res_array_U2 = cat(2,res_array_U2,sum(my_sim2.results.species_results_dict(i).value,2));
                        end
                end
                res_array_N_sum = sum(res_array_N,2);
                res_array_D_sum = sum(res_array_U,2);
                res_array_N2_sum = sum(res_array_N2,2);
                res_array_D2_sum = sum(res_array_U2,2);
                diff_curves = (res_array_N_sum+res_array_D_sum)-(res_array_N2_sum+res_array_D2_sum);
                mean_diff_curves = mean(diff_curves);
                capacity_krag = trapz(obj.results.T,diff_curves);
                capacity_krag_annual = capacity_krag/obj.results.T(end);
                fprintf('The krag capacity is: %8.0f.\n',capacity_krag);
                fprintf('The annual krag capacity is: %8.0f.\n',capacity_krag_annual);
                
                disp("Producing visual comparing launches vs no future launch case.")
                figure()
                plot(obj.results.T, sum(sum(obj.results.species_array(:,:,:),3),2), LineWidth=2.0, ...
                    DisplayName='With Launches')
                hold on
                plot(obj.results.T, sum(sum(my_sim2.results.species_array(:,:,:),3),2), LineWidth=2.0, ...
                    DisplayName='No Future Launches')
                xlabel('Time (years)')
                ylabel('Number of Objects')
                xlim([min(obj.results.T) max(obj.results.T)])
                title('Total Population')
                legend('Location','bestoutside', 'Interpreter', 'none')
                grid on

                disp("Producing visuals for Number-Time Product metric for Capacity.")
                figure
                hold on
                grid on
%                 set(gca,'FontSize',14)
                plot(obj.results.T,res_array_N2_sum+res_array_D2_sum,'b--','LineWidth',2)
                plot(obj.results.T,res_array_N_sum+res_array_D_sum,'m--','LineWidth',2)
                plot(obj.results.T,mean_diff_curves*ones(length(obj.results.T),1),'r-','LineWidth',2)
                plot(obj.results.T,diff_curves,'k--','LineWidth',2)
                area(obj.results.T,diff_curves,'FaceAlpha',0.4)
                xlabel('Time (years)')
                ylabel('Fragments')
                xlim([0 max(obj.results.T)])
                legend('No launches','With Launches','Mean','Difference','Location','bestoutside')
                title('Fragment-Year Index (Krag Capacity)')
    
                figure
                hold on
                grid on
%                 set(gca,'FontSize',14)
                plot(obj.results.T,diff_curves,'k--','LineWidth',2)
                plot(obj.results.T,mean_diff_curves*ones(length(obj.results.T),1),'r-','LineWidth',2)
                area(obj.results.T,diff_curves,'FaceAlpha',0.4)
                xlim([0 max(obj.results.T)])
                xlabel('Time (years)')
                ylabel('Fragments')
                legend('Difference (No Launch vs Launch)', 'Mean','Location','bestoutside')
                title('Fragment-Year Index (Krag Capacity)')
            else
                error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end
        end

        function cum_CSI = cum_CSI(obj)
            % Computes and displays cumulative Criticality of Space Index - CSI
            k = 0.6;
            life = @(h) exp(14.18*h^(0.1831)-42.94); % lifetime from power law fitting
            
            M_ref = 10000; % [kg]
            h_ref = 1000; % [km]
            life_h_ref = 1468; % [years] it corresonds to life0 = life(1000);
            
            D_ref = max(sum((obj.results.species_array(1,:,:)),3)./obj.scen_properties.V);

            den = M_ref*D_ref*life_h_ref*(1+k);
            
            cos_i_av = 2/pi; %average value of cosine of inclination in the range -pi/2 pi/2 calculated using integral average
            Gamma_av = (1-cos_i_av)/2;

            rgb_c = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],...
                [0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840]};
            if isprop(obj,"results")
                disp("Producing two visuals of CSI.")
                figure()
                hold on
                grid on
                CSI_S_sum_array = zeros(length(obj.results.T),0);
                CSI_D_sum_array = zeros(length(obj.results.T),0);
                for i2 = 1:length(obj.results.species_results_dict)
                    if i2>length(rgb_c)
                        colorset = [rand rand rand];
                    else
                        colorset = rgb_c{i2};
                    end
                    CSI_X_mat = zeros(length(obj.results.T),obj.scen_properties.N_shell);
                    name = obj.results.species_results_dict(i2).key;
                    if contains(name,'S') || contains(name,'D')
                        for i = 1:obj.scen_properties.N_shell
                            life_i = life((obj.scen_properties.R02(i)+obj.scen_properties.R02(i+1))/2);
                            num = life_i * (1+k*Gamma_av);
                            mass = obj.species_list(i2).species_properties.mass;
                            dum_X = mass * num;
                            D_X = obj.results.species_results_dict(i2).value(:,i)/obj.scen_properties.V(i2);
                            CSI_X_mat(:,i) = D_X*dum_X;
                        end
                        CSI_X_mat = CSI_X_mat./den;
                        CSI_X = sum(CSI_X_mat,2);
                        plot(obj.results.T,CSI_X,'DisplayName','CSI for '+strrep(name,'p','.'),'LineWidth',2,'Color',colorset)
                        if contains(name,'S') && ~contains(name,'D')
                            CSI_S_sum_array = cat(2,CSI_S_sum_array,CSI_X);
                        elseif contains(name,'D')
                            CSI_D_sum_array = cat(2,CSI_D_sum_array,CSI_X);
                        end
                    end
                end
            CSI_S_sum = sum(CSI_S_sum_array,2);
            CSI_D_sum = sum(CSI_D_sum_array,2);
%           threshold = max(CSI_S_sum+CSI_D_sum);
%           time_CSI = obj.results.T(1:end-1)';
%           available_capacity_S  = (threshold-CSI_S_sum(1:end-1,:))./(obj.results.T(end)-time_CSI);
%           available_capacity_D = (threshold-CSI_D_sum(1:end-1))./(obj.results.T(end)-time_CSI);
            plot(obj.results.T,sum(CSI_S_sum+CSI_D_sum,2),'DisplayName','Total CSI','LineWidth',2);
            xlabel('Time (years)')
            ylabel('CSI')
            title('CSI per Species Type')
            xlim([0 max(obj.results.T)])
            legend('Location','bestoutside', 'Interpreter', 'none')
            
            figure()
            grid on
            hold on
            plot(obj.results.T,CSI_S_sum,'DisplayName','Total CSI for Active Satellites','LineWidth',2);
            plot(obj.results.T,CSI_D_sum,'DisplayName','Total CSI for Derelict Satellites','LineWidth',2);
            plot(obj.results.T,sum(CSI_S_sum+CSI_D_sum,2),'DisplayName','Total CSI','LineWidth',2);
            xlabel('Time (years)')
            ylabel('Cumulative CSI')
            xlim([0 max(obj.results.T)])
            title('CSI for Active and Derelict Species')
            legend('Location','bestoutside', 'Interpreter', 'none')
            else
              error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for CSI computation.")
            end 
        end
        
        function init_vs_final = init_vs_final_vis(obj)
            %Produces a Figure comparing the total initial number of
            %objects and the total final number of objects 
            if isprop(obj,"results")
                species_array = obj.results.species_array;
                disp("Producing visual for the final vs initial populations.")
                alts = obj.scen_properties.HMid;
                figure()
                plot(alts, sum(species_array(1,:,:),3), LineWidth=2.0, ...
                    DisplayName='Initial Population')
                hold on
                plot(alts, sum(species_array(end,:,:),3), LineWidth=2.0, ...
                    DisplayName='Final Population')
                xlabel('Altitude (km)')
                ylabel('Number of Objects')
                xlim([min(alts) max(alts)])
                title('Total Population Initial vs Final')
                legend('Location','bestoutside', 'Interpreter', 'none')
                grid on
                set(gca, 'Yscale', 'log')
            else
              error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end 
        end

        function init_and_final = init_and_final_vis(obj)
            %Produces a Figure comparing the initial population of each 
            % species alongside a plot of the final population of each species 
            if isprop(obj,"results")
                disp("Producing visual for the initial and final species populations.")
                alts = obj.scen_properties.HMid;
                figure()
                for spec_ind = 1:length(obj.results.species_results_dict)
                    hold on
                    subplot(2,1,1)
                    plot(alts, obj.results.species_results_dict(spec_ind).value(1,:), LineWidth=2.0, ...
                    DisplayName=strrep(obj.results.species_results_dict(spec_ind).key,'p','.'))
                end
                xlabel('Altitude (km)')
                ylabel('Number of Objects')
                xlim([min(alts) max(alts)])
                title('Species Population Initial')
                legend('Location','bestoutside', 'Interpreter', 'none')
                grid on
                for spec_ind = 1:length(obj.results.species_results_dict)
                    hold on
                    subplot(2,1,2)
                    plot(alts, obj.results.species_results_dict(spec_ind).value(end,:), LineWidth=2.0, ...
                    DisplayName=strrep(obj.results.species_results_dict(spec_ind).key,'p','.'))
                end
                xlabel('Altitude (km)')
                ylabel('Number of Objects')
                xlim([min(alts) max(alts)])
                title('Species Population Final After Time = '+string(obj.results.T(end))+' Years')
                legend('Location','bestoutside', 'Interpreter', 'none')
                grid on
            else
              error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end 
        end

        function final_debris = final_debris_vis(obj)
            %Produces a Figure showing number of debris across altitudes 
            %by debris type 
            if isprop(obj,"results")
                disp("Producing visual for the final debris populations per type of debris.")
                figure()
                alts = obj.scen_properties.HMid;
                for i = 1:length(obj.results.species_results_dict)
                    name = obj.results.species_results_dict(i).key;
                    if contains(name,'N')
                        res_array = obj.results.species_results_dict(i).value;
                        plot(alts, res_array(end,:), LineWidth=2.0, ...
                        DisplayName=strrep(name,'p','.'))
                    end
                    hold on
                end
                xlabel('Altitude (km)')
                ylabel('Number of Debris')
                xlim([min(alts) max(alts)])
                title('Final Population of Each Debris Type')
                legend('Interpreter', 'none','Location','bestoutside')
                grid on
                set(gca,'YScale', 'log')
            else
              error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end 
        end

        function per_shell_indicator(obj, indicator_name, varargin)
            % Displays indicator values for indicators that have value per
            % shell.
            % Diff will take approximate numerical derivative of data.
            % units = (string, default ""), will be used to label units on graphs. Default is "" which displays no units.
            % take_deriv = (bool, default false), will approximate numerical derivative of
            % chosen indicator.
            % constraint = (double, default NaN), will use this value to define diverging 
            % color map for graph, centered on constraint value. Defaults
            % to NaN, which uses default colormap.
            % constraint_type = (string, default 'lower'). Valid
            % values are "lower" and "upper". If "lower", blue end of color
            % scale will be upper values.  If "upper", blue end of color
            % scale will be lower values.
            p = inputParser;
            addRequired(p, 'obj', @(x)class(x) == "simulation_class")
            addRequired(p, 'indicator_name')
            addOptional(p, 'units', '')
            addOptional(p, 'take_deriv', false, @islogical)
            addOptional(p, 'constraint', NaN)
            addOptional(p, 'constraint_type', 'upper', @(str) strcmpi(str,"upper") || strcmpi(str,"lower"))
            parse(p,obj, indicator_name, varargin{:});
            units = p.Results.units;
            take_deriv = p.Results.take_deriv;
            constraint  = p.Results.constraint;
            constraint_type = p.Results.constraint_type;

            % TODO: test diff, and test overlap between names for
            % integrated and non-integrated indicators.
            if isprop(obj,"results")
                alts = obj.scen_properties.HMid;
                disp("Producing visuals for the evolution of indicator.")
                figure()
                z = [];
                if take_deriv == false
                    if isfield(obj.results.indicators, indicator_name) && isfield(obj.results.integrated_indicators, indicator_name)
                        warning("Indicator ''" + indicator_name + "'' is present in both integrated indicators and indicators. Defaulting to indicator for display. Please rename one of the indicators to eliminate this warning.")
                    end
                    if isfield(obj.results.indicators, indicator_name)
                        z = obj.results.indicators.(indicator_name);
                        surf(obj.results.T,alts,z,'edgecolor','none')
                    elseif isfield(indicator_name, obj.results.integrated_indicators)
                        z = obj.results.integrated_indicators.(indicator_name);
                        surf(obj.results.T,alts,z.','edgecolor','none')
                    else
                        error("Indicator " + indicator_name + " was not found. Please check the name and try again")
                    end
                elseif take_deriv == true
                    if isfield(obj.results.indicators, indicator_name) && isfield(obj.results.integrated_indicators, indicator_name)
                        warning("Indicator ''" + indicator_name + "'' is present in both integrated indicators and indicators. Defaulting to indicator for display. Please rename one of the indicators to eliminate this warning.")
                    end
                    if isfield(obj.results.indicators, indicator_name)
                        approx_diff = diff(obj.results.indicators.(indicator_name), 1, 2) ./diff(obj.results.T, 1, 1).';
                        Tdiff = (obj.results.T(1:end-1,:) + obj.results.T(2:end,:)) / 2;
                        z = approx_diff;
                        surf(Tdiff,alts,z,'edgecolor','none')
                    elseif isfield(obj.results.integrated_indicators, indicator_name)
                        approx_diff = diff(obj.results.integrated_indicators.(indicator_name), 1, 2) ./diff(obj.results.T, 1, 1).';
                        Tdiff = (obj.results.T(1:end-1,:) + obj.results.T(2:end,:)) / 2;
                        z = approx_diff;
                        surf(obj.results.T,alts,z,'edgecolor','none')
                    else
                        error("Indicator " + indicator_name + " was not found. Please check the name and try again")
                    end
                end

                if ~isnan(constraint)
                    if constraint_type == "upper"
                        cmap = flipud(lbmap(256,'BrownBlue'));
                    elseif constraint_type == "lower"
                        cmap = lbmap(256,'BrownBlue');
                    end
                    zdiff_max = abs(constraint - max(z, [], "all"));
                    zdiff_min = abs(constraint - min(z, [], "all"));
                    clim_bounds = [constraint-max(zdiff_min, zdiff_max), constraint + max(zdiff_min, zdiff_max)];
                    colormap(cmap)
                    %colormap(centered)
                    c = colorbar();
                    clim(clim_bounds)
                else
                    %colormap(hot)
                    c = colorbar();
                end

                c.Label.Interpreter = 'none';
                if units == ""
                    c.Label.String = indicator_name;
                    zlabel(indicator_name, 'Interpreter', 'none')
                else
                    c.Label.String = indicator_name + " ( " + units + " )" ;
                    zlabel(indicator_name + " ( " + units + " )", 'Interpreter', 'none')
                end
                    
                xlabel('Time (years)', "FontSize", 14)
                ylabel('Altitude (km)', "FontSize", 14)
                xlim([0 max(obj.results.T)])
                ylim([min(alts) max(alts)])
                title(indicator_name, 'Interpreter', 'none')
                grid on
            else
                error("Simulation does not contain results. Please run the function run_model(obj, x0) to produce simulation results required for visualisations.")
            end
        end

    end %end methods
end %end class

