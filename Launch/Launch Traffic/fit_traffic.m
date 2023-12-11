function [Lambda_Func] = fit_traffic(times, launches, launch_period, graph)
    % Construct a function composed of Gaussians to create a function from
    % discrete launch rates.
    % times: a 1 by n vector of doubles for time in years
    % launches: a 1 by n vector of satellites launched at the times in the
    % times vector
    % launch_period: the period of time that 4-sigma (99.994%) of launches
    % should occur. Default is 1 (year).
    % graph: bool. Displays graphs if true.

    if ~isempty(launch_period)
        launch_period = 1;
    end

    if ~isempty(graph)
        launch_period = 1;
    end


    x = launches;
    tspan0 = times;
    models = {};
    
    % If no launches in a window, make curve 0.
    if isempty(launches) == 1
        Lambda_Func = @(t) 0*t;
        return
    end
    for xind = 1:length(x)
        peak = x(xind);
        mu = tspan0(xind);
        %sigma = sqrt(2)/(2*peak*sqrt(pi)); %So that the max value of the normal will be lambdasigma = sqrt(2)/(2*peak*sqrt(pi)); %So that the max value of the normal will be lambda
        
%         if xind == 1 % First point
%             %disp("First Point")
%             x_range = linspace(tspan0(xind) - launch_period/2, tspan0(xind+1) + launch_period/2, 100);
%             %x_range = linspace(tspan0(xind), tspan0(xind+1), 100);
%         elseif xind == length(x)
%             %disp("Last Point")
%             x_range = linspace(tspan0(xind-1), tspan0(xind), 100);
%         else
%             x_range = linspace(tspan0(xind-1), tspan0(xind+1), 100);
%         end
        %scale_factor = 1;
        
        % Scale the normal distribution so 4 sigma within the launch_period
        sigma = launch_period/8;
        fun = @(scale_factor) (scale_factor * integral(@(x)normpdf(x,mu,sigma), min(tspan0), max(tspan0)) - peak);
        [scale_factor] = fzero(fun,[0.0 1e6]);
        %[scale_factor] = fzero(fun,[0.0 1e6],optimset("Display","iter"));
        
        %fun = @(x) -sqrt(pi)*scale_factor*erf(-peak*sqrt(pi)*x/scale_factor + mu*peak*sqrt(pi)/scale_factor)/(2*sqrt(pi)) - scale_factor;
        %scale_factor = @(x) integral(@(x) normpdf(x,mu,sigma), min(tspan0), max(tspan0));
        %scale_factor = x(xind)/integral(@(x) normpdf(x,mu,sigma), min(tspan0), max(tspan0));
        
        models{xind} = @(x) scale_factor * normpdf(x,mu,sigma);
        %Lambda_Func = @(x) sum(reshape(cell2mat(arrayfun(@(idx) models{idx}(x), 1:xind, 'uniform', 0)), xind, []));
        %models = sum(models)(tspan(1:xind));
        %disp(["xind:", xind, "mu:", mu, "sigma:", sigma, "lambda", x(xind), "int", integral(models{xind}, 0, 100), "cumsum", x(xind)])
        %disp(["tspan", tspan0])
        %disp(["x_range", x_range])
        
        % Individual graphs
%         if graph
%             %Individual Fit Functions To Now
%             figure(1)
%             clf
%             hold on
%             plot(tspan0(1:xind),x(1:xind), 'ok', "DisplayName", "Data")
%             f_xint = linspace(0, max(tspan0(1:xind) + launch_period/2), 100);
%             for pind = 1:length(models)
%                 plot(f_xint, models{pind}(f_xint), ":", 'HandleVisibility','off')
%             end
%             legend()
%             hold off
%             
%             % Only Latest
%             launch_start = tspan0(xind) - launch_period;
%             launch_end =  tspan0(xind) + launch_period;
%             x_range = linspace(launch_start, launch_end, 100);
%             fig = figure(2);
%             clf
%             hold on
%             plot(tspan0(xind),x(xind),'ok', "DisplayName", "Data");
%             plot(x_range, models{xind}(x_range), ":", "DisplayName", "Model");
%             hold off
%         end
    end % end loop over points
  
    Lambda_Func = @(t)sum(reshape(cell2mat(arrayfun(@(idx) models{idx}(t), 1:length(x), 'uniform', 0)).', length(t),[]).');

    if graph
        % Plot of launch curves and aggregate function
        f_full = linspace(0, max(tspan0) + launch_period, 200);
        fig = figure(4);
        clf
        hold on
        %sums = zeros(size(f_full));
        for pind = 1:length(models)
            plot(f_full, models{pind}(f_full), ":",'HandleVisibility','off')
            %sums = sums + models{pind}(f_full);
        end
        plot(tspan0,x, 'or', "DisplayName", "Fit Function") %Data
        % plot(f_full,sums, "-r", "DisplayName", "Sums")
        plot(f_full,Lambda_Func(f_full), "--b", "DisplayName", "Fit Function")
        legend
        xlabel('Time (years)')
        ylabel('Number of Satellites Launched')
        hold off

        % Compare the cumulative sum and integral
        fig = figure(5);
        clf
        hold on
        cumInt = @(t)integral(Lambda_Func, 0, t);
        plot(tspan0,cumsum(x), "--b", "DisplayName", "Cumulative Sum")
        plot(f_full,arrayfun(cumInt,f_full), "-r", "DisplayName", "Fit Func. Integral")
        hold off
        legend('Location',"northeast")
        xlabel('Time (years)')
        ylabel('Cumulative Number of Satellites Launched')
        
    end
end