function [Fdot] = drag_func(t, species_properties, scen_properties)

    Fdot = zeros(scen_properties.N_shell, 1, 'sym');

    if species_properties.drag_effected == true
        h = scen_properties.R02;
        rho = scen_properties.dens_model(t, h, scen_properties);
        
        for k=1:scen_properties.N_shell
            if k<scen_properties.N_shell % Not top shell
                n0=(species_properties.sym(k+1));
                h = scen_properties.R02(k+1);
                rho_k1 = rho(k+1);
                rvel_upper = -rho_k1*species_properties.beta*sqrt(scen_properties.mu*scen_properties.R0(k+1))*(24*3600*365.25);% drag flux
            else
                % ASSUMPTION: No flux coming down from highest shell.
                n0 = 0;
                h = scen_properties.R02(k+1);
                rho_k1 = rho(k+1);
                rvel_upper =-rho_k1*species_properties.beta*sqrt(scen_properties.mu*scen_properties.R0(k+1))*(24*3600*365.25);% drag flux
            end % End k<scen_properties.N_shell
            
            rhok = rho(k);
            rvel_current=-rhok*species_properties.beta*sqrt(scen_properties.mu*(scen_properties.R0(k)))*(24*3600*365.25);
            Fdot(k, 1) = +n0*rvel_upper/scen_properties.Dhu + rvel_current/scen_properties.Dhl * species_properties.sym(k);
            
        end % End k=1:scen_properties.N_shell
    else % Not drag_effected
        for k=1:scen_properties.N_shell
            Fdot(k,1) = 0 * species_properties.sym(k);
        end % End k=1:scen_properties.N_shell
    end % End species_properties.drag_effected == true

    Fdot = sym2cell(Fdot);


    
    