function outmatsats = mega2matsats(constellationFile,constellationSheet, scenarioYrs, time0)
% from launch xlsx file (eg megaconstellationLaunches.xlsx) to an 
% instantiation of matsats for X yrs

% Example: 
% constellationFile = 'megaconstellationLaunches.xlsx'; 
% constellationSheet = 'forexport2';  % fewer constellations

getidx;

% m = readmatrix(fname);

% from setupTLE_LNT:
s = detectImportOptions(constellationFile,'Sheet',constellationSheet);
tab = readtable(constellationFile,s);

% fill empty entries
for ind = 1:size(tab,1)
    if isnan(tab(ind,:).FirstLaunch)
        tab(ind,:).FirstLaunch = 2030;  % if no launch specified, start launch in 2030
        tab(ind,:).FinishLaunch = 2045;  % if no launch specified, end launch in 2045
    end
    if isnan(tab(ind,:).FinishLaunch)        
        tab(ind,:).FinishLaunch = tab(ind,:).FirstLaunch  + 10;  % if no end is specified, end launch in 15 years
    end
    if isnan(tab(ind,:).mass)
        tab(ind,:).mass = tab(1,:).mass;  % if no mass, copy Starlink-v1's        
    end
    if isnan(tab(ind,:).radius)
        tab(ind,:).radius = (tab(ind,:).mass / tab(1,:).mass)^(1/3) * tab(1,:).radius;  % use Starlink-v1's density
    end
end

tab.missionlife = ones(size(tab,1),1) * 8;      % add missionlife as 8 yrs


% from main_mc_LNT:
constellation = tab;
% for ind = 1:size(constellation,1)
%     if year(current_time) >= constellation.FirstLaunch(ind) && year(current_time) <= constellation.FinishLaunch(ind)
%         % if still in ramp-up phase
%         Nrampperyear(ind) = (constellation.TotalSatsPlanned(ind) - constellation.SatsOnStn(ind)) / (constellation.FinishLaunch(ind) - year(time0) + 1);
%         % maintenance launches needed: maintain by launching Num-in-orbit / missionLife [not exact]
%         Nmaintperyear(ind) = Nrampperyear(ind) * (year(current_time) - constellation.FirstLaunch(ind)) / constellation.missionlife(ind);  
%     elseif  year(current_time) >= constellation.FirstLaunch(ind)
%         % if constellation is done filling up
%         Nrampperyear(ind) = 0;
%         Nmaintperyear(ind) = constellation.TotalSatsPlanned(ind) / constellation.missionlife(ind);
%     else
%         Nrampperyear(ind) = 0;
%         Nmaintperyear(ind) = 0;
%     end
% end


% loop over mc
YEAR2DAY = 365.2425; 
dt_days = 5; 
earthRadius = 6371000; % m
outmatsats = [];
constLaunch = [];

for d = 1:dt_days:scenarioYrs*YEAR2DAY
    current_time = time0 + d;
    jd = juliandate(current_time);

    for ind = 1:size(constellation,1)
        if year(current_time) >= constellation.FirstLaunch(ind) && year(current_time) <= constellation.FinishLaunch(ind)
            % if still in ramp-up phase
            Nrampperyear(ind) = (constellation.TotalSatsPlanned(ind) - constellation.SatsOnStn(ind)) / (constellation.FinishLaunch(ind) - year(time0) + 1);
            % maintenance launches needed: maintain by launching Num-in-orbit / missionLife [not exact]
            Nmaintperyear(ind) = Nrampperyear(ind) * (year(current_time) - constellation.FirstLaunch(ind)) / constellation.missionlife(ind);
        elseif  year(current_time) >= constellation.FirstLaunch(ind)
            % if constellation is done filling up
            Nrampperyear(ind) = 0;
            Nmaintperyear(ind) = constellation.TotalSatsPlanned(ind) / constellation.missionlife(ind);
        else
            Nrampperyear(ind) = 0;
            Nmaintperyear(ind) = 0;
        end
    end



    if size(constellation,1) > 0
        % debug
        %disp(d)
        if mod(d, dt_days*100+1) == 0
            fprintf('\t\t %i: Number of constellations ramping: %i \t\t maintaining: %i \n', ...
                year(current_time), ...
                sum( year(current_time) >= constellation.FirstLaunch & year(current_time) <= constellation.FinishLaunch), ...
                sum( year(current_time) > constellation.FinishLaunch));
        end
        avglaunch = (Nmaintperyear + Nrampperyear)/(YEAR2DAY/dt_days);  % avg launch per timestep (fractional)
        intlaunch = floor(avglaunch);  % guaranteed launch amount per time step (int)
        addlaunch = rand(size(avglaunch)) < (avglaunch - intlaunch);    % additional launches - sampled from fractional part
        curConstLaunch = intlaunch + addlaunch;
    
        constMatsat = nan(size(avglaunch,2), 24);
        constMatsat(:,[idx_a, idx_ecco, idx_inclo, idx_constel, idx_bstar, idx_mass, idx_radius, idx_missionlife]) = ...
            [constellation.Altitude / earthRadius*1000 + 1, ones(size(avglaunch))'/1e6, ...
            deg2rad(constellation.Inclination), ones(size(avglaunch))', ...
            0.5 * 2.2 * constellation.radius .^2 ./ constellation.mass * 0.157, ...
            constellation.mass, constellation.radius, constellation.missionlife];
        constMatsat(:, idx_a_desired) = constMatsat(:, idx_a);
        constMatsat(:,[idx_controlled, idx_objectclass]) = [ ones(size(avglaunch))',  ones(size(avglaunch))'];
    
        constLaunch = constMatsat(repelem(1:size(constMatsat,1), curConstLaunch), :  ); % populate constellation launch matsats
    
        % bstar = 0.5 * 2.2 * ipop(:, idx_radius).^2 ./ ipop(:, idx_mass) * 0.157;
    
        % Scramble launch time within this time step (dt_days), argpo, mo, nodeo
        constLaunch(:,idx_launch_date) = jd + rand(size(constLaunch,1),1) * dt_days;
        constLaunch(:,[idx_argpo, idx_mo, idx_nodeo]) = 2 * pi * rand(size(constLaunch,1),3);
    
        outmatsats = [outmatsats; constLaunch];
%         param.maxID = param.maxID+size(out_future,1)+size(constLaunch,1);       %update maximum ID in population
%         count_tot_launches = count_tot_launches+size(out_future,1)+size(constLaunch,1);
%         out_future = [out_future; constLaunch];
    end
end


end