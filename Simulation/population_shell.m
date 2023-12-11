function [xdot] = population_shell(t, x, obj)
    % For time-varying atmosphere, density needs to be computed within the
    % integrated function, not as an argument outside it.
    % Inputs:
    %   t is a time in years from start date.
    %   X is the equation state
    %   obj is a simulation object.
    %   var is the list if state variables assembled in run_model
    
    obj.scen_properties.X = x;
    obj.scen_properties.t = t;
    %disp("X")
    %disp(x.')
    
    for fun_i = 1:length(obj.scen_properties.functions_to_run_each_model_loop)
        %disp(fun_i)
        fun = obj.scen_properties.functions_to_run_each_model_loop{fun_i};
        fun(x, t);
    end
    %cellfun(@(fun) fun(x, t),functions_to_run_each_model_loop,'UniformOutput',false);

    if obj.scen_properties.disp_times
        disp("t:" + num2str(t))
    end
    
    % %TODO: pull this out for when static expotential model. Will be faster.
    xdot = obj.xdot_fun(x);
    %disp(t)

    %% Launch rates
    % Use non-symbolic if possible since it is much faster.
    if obj.scen_properties.sym_lambda == false
        full_lambda = obj.scen_properties.full_lambda;
        full_lambda_val = cellfun(@(fun) fun(x, t),full_lambda,'UniformOutput',false);
        try
            %full_lambda_val = double([full_lambda_val{:}]).'; %cast any 0 symbolics to double 0
            %xdot = xdot + full_lambda_val;
            xdot = xdot + [full_lambda_val{:}].'; %lambda_fun(x.');
            %for i = 1:length(full_lambda_val)
            %    disp(class(full_lambda_val{i,1}))
            %end
    
        catch
            %TODO: fallback to slower symbolic.
            lambda_fun = matlabFunction(full_lambda_val,'Vars',obj.var_col);
            xdot = xdot + lambda_fun(x.');
        end
        
        %xdot = obj.xdot_fun(x);
    end
    
    if obj.scen_properties.time_dep_density == true
        rho = obj.scen_properties.dens_model(t, obj.scen_properties.R02, obj.scen_properties);
        
        rho_mat_k = repmat(rho(1:end-1), width(obj.scen_properties.species), 1);
        rho_mat_kp1 = repmat(rho(2:end), width(obj.scen_properties.species), 1);
        
        if isprop(obj.scen_properties, 'indicator_var_list')
            rho_mat_k = [rho_mat_k; zeros(obj.scen_properties.num_integrated_indicator_vars, 1)]; %Pad for indicator variables
            rho_mat_kp1 = [rho_mat_kp1; zeros(obj.scen_properties.num_integrated_indicator_vars, 1)]; %Pad for indicator variables
        end
        
        full_drag = obj.drag_term_upper(x) .* rho_mat_kp1 + obj.drag_term_cur(x) .* rho_mat_k;
        xdot = xdot + full_drag;
        %sprintf('%0.16f',t)
    else
        xdot = xdot;
        %sprintf('%0.16f',t)
    end

%     if obj.scen_properties.time_dep_density == true
%         full_drag = zeros(obj.scen_properties.N_shell, width(obj.species_list), 'sym');
%         for i = 1:length(obj.species_list)
%             Fdot = obj.species_list(i).drag_func(t, obj.species_list(i).species_properties, obj.scen_properties);
%             full_drag(:,i) = Fdot;
%         end
%     
%         full_drag_col = reshape(full_drag, [obj.scen_properties.N_shell*width(obj.scen_properties.species),1]);
%         %full_drag_f = matlabFunction(full_drag_col,'Vars',obj.var_col);
%         full_drag_f = double(subs(full_drag_col, obj.var_col, x));
%         xdot = obj.xdot_fun(x); % + full_drag_f;
%         sprintf('%0.16f',t)
%     else
%         xdot = obj.xdot_fun(x);
%         
%         sprintf('%0.16f',t)
%     end

end