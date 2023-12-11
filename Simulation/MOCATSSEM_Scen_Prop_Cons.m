function [scen_properties] = MOCATSSEM_Scen_Prop_Cons(varargin)
    % This function assembers the scen_properties object with some basic
    % default values
    % Inputs: fields and their purposes should be self-documenting.
    % Outputs: scen_properties object.
    
    p = inputParser;
    addOptional(p,'start_date', datetime(2022,1,1,0,0,0), @isdatetime);
    addOptional(p,'simulation_duration', 200, @isfloat);
    addOptional(p,'steps', -1, @isfloat); % default to simulation duration
    addOptional(p,'density_model', @static_exp_dens_func, @(x) isa(x,'function_handle') );
    addOptional(p, 'integrator', @ode15s, @(x) isa(x,'function_handle'))
    addOptional(p, 'shells', 20, @isint)
    addOptional(p, 'max_altitude', 1400, @isint)
    addOptional(p, 'min_altitude', 200, @isint)
    addOptional(p, 'density_filepath',[],@isstring)
    addOptional(p, 'delta',10,@isint)
    addOptional(p, 'LC', .1,@isfloat)
    addOptional(p, 'v_imp',10,@isfloat)
    
    parse(p, varargin{:})
    VAR.start_date = p.Results.start_date;
    VAR.start_year = year(VAR.start_date);
    VAR.start_month = month(VAR.start_date);
    VAR. start_day = day(VAR.start_date);
    VAR.start_hour = hour(VAR.start_date);
    VAR.start_minute = minute(VAR.start_date);
    VAR.start_second = floor(second(VAR.start_date));
    VAR.start_second_fraction = second(VAR.start_date) - floor(second(VAR.start_date));
    
    VAR.simulation_duration = p.Results.simulation_duration;
    VAR.N_steps = p.Results.steps;
    if VAR.N_steps == -1
       VAR.N_steps = round(VAR.simulation_duration);
    end
    
    VAR.N_shell = p.Results.shells;
    VAR.h_max = p.Results.max_altitude;
    VAR.h_min = p.Results.min_altitude;
    VAR.LC = p.Results.LC; % minimum size of fragments [m]
    VAR.v_imp = p.Results.v_imp;
    VAR.delta = p.Results.delta;
    
    VAR.dens_model = p.Results.density_model;
    if functions(VAR.dens_model).function == "static_exp_dens_func"
        VAR.time_dep_density = false;
    elseif functions(VAR.dens_model).function == "JB2008_dens_func"
        VAR.time_dep_density = true;
        if isempty(p.Results.density_filepath)
            VAR.density_filepath = fullfile(".","Atmosphere Model", "JB2008", "Precomputed", "dens_highvar_2000.mat");
        end
    else
        warning("Unable to detect whether atmosphere model " + functions(VAR.dens_model).function + ...
                " is time dependent. Please manually set " + ...
                "scen_properties.time_dep_density prior to building model.")
    end
    
    VAR.integrator = p.Results.integrator;
    variable_step_integrators = ["ode45", "ode21", "ode113", "ode78", "ode89", ...
                                 "ode15s", "oder23s", "ode23t", "ode23tb"];
    fixed_step_integrators = ["ode1", "ode2", "ode3", "ode4", "ode5"];
    integrator_name = functions(VAR.integrator).function;
    if ismember(integrator_name, variable_step_integrators)
        VAR.integrator_type = 'variable';
    elseif ismember(integrator_name, fixed_step_integrators)
        VAR.integrator_type = 'fixed';
    else
        warning("Unable to detect whether integrator " + integrator_name + ...
                " is fixed or variable step. Please manually set " + ...
                "scen_properties.integrator_type prior to building model.")
    end
    
    
    % Establish parameter-dependent quantities and constants
    VAR.scen_times = linspace(0,VAR.simulation_duration,VAR.N_steps)';
    
    % Constants
    rad = pi/180;
    mu_km=398600.4415;%km^3/s^2
    mu=3.986004418e14;%meters^3/s^2
    VAR.re = 6378.1366; % [km]
    years_s = 365*24*3600;
    
    % Parameters needed for all MOCAT-SSEMs
    %VAR.N_step = 5; % how many steps to calculate between t0 and tf
    R0=linspace(VAR.h_min,VAR.h_max,VAR.N_shell+1);
    R02=linspace(VAR.h_min,VAR.h_max,VAR.N_shell+1);
    VAR.HMid = R02(1:end-1) + diff(R02)/2;
    VAR.deltaH=R02(2)-R02(1); % thickness of the shell [km]
    
    %VAR.v_imp = 10; % impact velocity [km/s]
    R1=R0;
    R0=(VAR.re+R0)*1000;
    % V=4*pi*R0.^2*deltaH*1000;
    VAR.V=4/3*pi*(diff(R0.^3)); % volume of the shells [m^3]
    VAR.v_imp2 = VAR.v_imp * ones(size(VAR.V)); % impact velocity [km/s] %Shell-wise
    
    VAR.v=VAR.v_imp2*1000*(24*3600*365.25);% impact velocity [m/year]
    %VAR.v2=VAR.v_imp*1000*(24*3600*365.25);% impact velocity [m/year]
    VAR.Dhl=VAR.deltaH*1000;
    VAR.Dhu=-VAR.deltaH*1000;
    
    %Integration stuff
    options = odeset('reltol', 1.e-4,'abstol', 1.e-4);
    VAR.options = options;
    VAR.mu = mu;
    VAR.R0 = R0;
    VAR.R02 = R02;
    
    scen_properties = scen_properties_class(VAR);

end

function [bool,idx] = isint(x)
    % from https://www.mathworks.com/matlabcentral/answers/16390-integer-check
    % Check whether input is integer or not
    % Inf and NaN are not integers
    if ~isnumeric(x)
       error('Input must be a numeric, not a %s.',class(x))
    end
    bool = (mod(x,1) == 0);
%     bool = (round(x) == floor(x));  % Other approach. Fails with Inf and
%     NaN.
    idx = find(bool); 
% Manolín Sept-2019
end
