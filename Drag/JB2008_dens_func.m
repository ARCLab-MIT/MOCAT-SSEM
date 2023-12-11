function [p] = JB2008_dens_func(t, h, scen_properties)
    % JB2008 density wrapper. Will error if not called for values of t and
    % h interpolated during build_model.
    %   t is time from scenario start in years (scalar only)
    %   h is the height above ellipsoid in km. (array [only if model altitudes] or scalar)
    %   scen_properties is a structure with properties for the scenario
    %       scen_properties.density_filepath is the path to the precomputed
    %       density file created by dens_prediction or included in the
    %       precomputed folder.
    %       scen_properties.interp_dens_var has the interpolated data for
    %       relevant t and h periods.
    %   density p is returned in kg/km^3.
    
    if ~isprop(scen_properties, 'density_filepath')
        error("scen_properties.density_filepath must be provided to " + ...
            "use JB2008 density information.")
    end

    if ~isprop(scen_properties, 'interp_dens_var')
        error("scen_properties.interp_dens_var must be computed in " + ...
              "simulation_class.build_model before densities can be called")
    end

    if class(t) == "sym"
        warning("Symbolic t passed to time-varying atmospheric model. Assuming constant density value based on t = 0.")
        t = 0;
    end

    % Find t index for scen_properties.interp_dens_var col
    closest = interp1(scen_properties.scen_times,scen_properties.scen_times,t,'nearest');
    t_idx = find(scen_properties.scen_times == closest);
    %disp("JB2008 t_idx: " + num2str(t_idx))
    
    if all(size(h) == size(scen_properties.R02)) && all(h == scen_properties.R02)
        % If you call with the same altitudes as the model, skip interp. 
        p = scen_properties.interp_dens_var(:, t_idx);
        return
    elseif strcmp(class(h), 'double')
        % If it's a fixed altitude, do the interpolation.
        h_idx = find(scen_properties.R02 == h);
        %Use closest neighbor if can't find, 0 if outside range.
        if isempty(t_idx) || isempty(h_idx) 
            if h < scen_properties.h_min
                error(num2str(h) + "is below minimum interpolated height," ...
                      + num2str(scen_properties.h_min) + ".")
            elseif h > max(scen_properties.R02)
                error(num2str(h) + "is above maximum interpolated height," ...
                      + num2str(max(scen_properties.R02)) + ".")
            elseif t < min(scen_properties.scen_times) || t > max(scen_properties.scen_times)
                error(num2str(t) + "is outside interpolated timespan" ...
                      + num2str([min(scen_properties.scen_times), ...
                      max(scen_properties.scen_times)]) + ".")
            end
            p = interp2(scen_properties.scen_times.', scen_properties.R02.', ...
                        scen_properties.interp_dens_var, t, h, 'nearest', 0);
        else
            p = scen_properties.interp_dens_var(h_idx, t_idx);
        end
    end


