function [upper_term, current_term] = drag_func_sym(species_properties, scen_properties)
    
    % This function makes the symbolic portions of the drag function
    % (without density). This allows for time-varying rhos to be used at
    % much better speed.
    sec_per_year = 24*3600*365.25;

    rvel_upper = zeros(scen_properties.N_shell, 1, 'sym');
    rvel_current = zeros(scen_properties.N_shell, 1, 'sym');
    upper_term = zeros(scen_properties.N_shell, 1, 'sym');
    current_term = zeros(scen_properties.N_shell, 1, 'sym');

    if species_properties.drag_effected == true      
        for k=1:scen_properties.N_shell
            %disp(k)
            if k<scen_properties.N_shell % Not top shell
                n0=(species_properties.sym(k+1));
                rvel_upper(k) = -species_properties.beta*sqrt(scen_properties.mu*scen_properties.R0(k+1))*sec_per_year;% drag flux
            else
                % ASSUMPTION: No flux coming down from highest shell.
                n0 = 0;
                rvel_upper(k) =-species_properties.beta*sqrt(scen_properties.mu*scen_properties.R0(k+1))*sec_per_year;% drag flux
            end % End k<scen_properties.N_shell
            
            rvel_current(k) =-species_properties.beta*sqrt(scen_properties.mu*(scen_properties.R0(k)))*sec_per_year;
            upper_term(k) = n0*rvel_upper(k)/scen_properties.Dhu;
            current_term(k) = rvel_current(k)/scen_properties.Dhl * species_properties.sym(k);
            
        end % End k=1:scen_properties.N_shell
    else % Not drag_effected
        % Keeps zeros from initialization.
    end % End species_properties.drag_effected == true
    %upper_term = sym2cell(upper_term);
    %current_term = sym2cell(current_term);
    


    
    