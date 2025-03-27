function cfgMC = setupMCeqDynamic(rngseed, nyears, ...
                    R02, N_shell, diameter, mass, xopt, alpha, alpha_active, P, Dt, Area)

x = reshape(xopt,N_shell,4); 

Seq = [round(x(:,1))];    % x = reshape(xopt,14,4); round(x(:,1))
Deq = [round(x(:,2))];    % round(x(:,2))
Neq = [round(x(:,3))];
lam = round(Seq ./ Dt);

iPop = round([Seq,Deq,Neq]);

% Set conversion units
cfgMC.DAY2MIN = 60*24;
cfgMC.DAY2SEC = cfgMC.DAY2MIN*60;
cfgMC.YEAR2DAY = 365.2425;  % days(years(1)) https://www.mathworks.com/help/matlab/ref/duration.years.html
cfgMC.YEAR2MIN = cfgMC.YEAR2DAY * cfgMC.DAY2MIN;
cfgMC.rad = pi/180;

% GLOBAL VARIABLES
whichconst = 84;
[tumin, mu_const, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc( whichconst );
cfgMC.tumin = tumin;
cfgMC.mu_const = mu_const;
cfgMC.radiusearthkm = radiusearthkm;
cfgMC.j2 = j2;
cfgMC.omega_earth = 2*pi/(cfgMC.DAY2SEC);

% SCENARIO PARAMETERS
cfgMC.PMD = P;                   % post mission disposal
cfgMC.alph = alpha;                     % active control failure probability with one sat active
 cfgMC.alph_a = alpha_active;                   % active control failure probability with both active
    cfgMC.alph_c_intra = 0;           % active control failure  between objects of different const
    cfgMC.alph_c_inter = 0;           % active control failure  between objects of same constellation
cfgMC.orbtol = 7.5;                 % orbital tolerance for controlled satellites [km]
cfgMC.step_control = 1;             % timesteps to check for orbit control tolerance
cfgMC.P_frag = 0;                   % probability per day of Rocket Body Fragmentation (if P_frag=0, explosions are not considered!)
cfgMC.P_frag_cutoff = inf;          % age at which objects don't explode.   ESA = 18 yrs.
cfgMC.altitude_limit_low = R02(1);     % lower limit of altitude [km]
cfgMC.altitude_limit_up = R02(end);     % upper limit of altitude [km] % IADC: perigee < 2000
cfgMC.missionlifetime = Dt;          % [years] operational life of payloads

% paramSSEM structure needed to summarize data into bins for TLE
% simulation analysis.  Needed for Fast_MC2SSEM_population.m
% Below snippet from prepare_pop_in.m
paramSSEM.N_shell = N_shell;
paramSSEM.h_min = cfgMC.altitude_limit_low;  % [km] defined in setup above
paramSSEM.h_max = cfgMC.altitude_limit_up;
% paramSSEM.R02 = linspace(paramSSEM.h_min,paramSSEM.h_max,paramSSEM.N_shell+1);
paramSSEM.R02 = R02;
paramSSEM.re = radiusearthkm;   % km

% PAD to fill shells not defined in iPop
    paramSSEM.lam2 = zeros(paramSSEM.N_shell, 1);
    paramSSEM.lam2(1:size(lam,1)) = lam;

    paramSSEM.iPop = zeros(paramSSEM.N_shell, size(iPop,2));
    paramSSEM.iPop(1:size(iPop,1),1:size(iPop,2)) = iPop;

paramSSEM.lam1 = paramSSEM.lam2(:,1);                % just the launched payload

paramSSEM.radius = diameter./2;     % m
paramSSEM.mass = mass;              % kg
paramSSEM.A = Area;                 % m^2

cfgMC.paramSSEM = paramSSEM;

% SET TIMES
cfgMC.time0 = datetime(2023,1,1);
t0_prop = 0;                            % [min]
nyears = nyears;                        % [years]
tf_prop = cfgMC.YEAR2MIN * nyears;      % [min]
cfgMC.dt_days = 5;                      % [days] sampling time for CUBE method and propagation
DeltaT = cfgMC.dt_days*cfgMC.DAY2MIN;     % [min] sampling time for CUBE method and propagation
cfgMC.tsince = t0_prop:DeltaT:t0_prop+tf_prop; % [min] set the length of the simulation, 200 years is typical;
cfgMC.n_time = length(cfgMC.tsince);

% LAUNCHES
Simulation = 'TLE';                     % 'TLE'
cfgMC.launch_model = 'matsat';
cfgMC.total_launch_per_year = 0;        % only for 'TLE' with TLElaunchRepeat = 0 for 'random' launch_model; used in launches_current_step_vec.m 'random' and future_traffic_model_vec;
cfgMC.launch_increase_per_year = 0;     % increase in launch rate per year since t0
TLElaunchRepeat = 0;                    % 0: 'random' launch via poisson distribution (see initSim below)
% 1: ESA style launch (repeat launches between years X and Y)
cfgMC.launchRepeatYrs = [0,1];    % Min/max year of obj to repeatedly launch (inclusive years)
                                % Only used if TLElaunchRepeat == 1
                                % Range of this used to repeat launch (!!)
                                % the date of repeatlaunches (repeatlaunches(:,idx_launchdate) 
                                % does not matter!! 
% ESA: 2001-2009 repeated (~ 75 / yr via Fig 2.14 in ESA 2022)
cfgMC.launchRepeatSmooth = 1;           % [0/1] average out the above so yearly launch rate remains the same
cfgMC.launchMinPerShell = 0;           % [integer >= 0] min number of launches (payload) per shell as specified in paramSSEM.R02
% cfgMC.xLaunch = xLaunch;                % multiplier for launched population

        % Constellation
        cfgMC.constellationFile = ''; 
        cfgMC.constellationSheet = '';
        cfgMC.constellation = ''; % = setupConstellation(cfgMC);  % populate cfgMC.constellation (table)


% Modify initial population
% Fill in missing satellite physical parameters
cfgMC.fillMassRadius = 1;       % 0: don't fill in missing DISCOS data (many objects with 0 radius and/or mass)
                                % 1: ESA's method -- assume spherical aluminum depending on RCS size (S/M/L)
                                % 2: resampling method
% cfgMC.xPop = xPop;              % multiplier for initial population
cfgMC.noRBs = 1;                % [0/1] remove RBs from simulation (helpful for MOCAT-3 validation)
cfgMC.physicalBstar = 1;        % [0/1] recalculate B* as Bstar = 1/2 * Cd * A/m * 0.157e6
                                % see line 89 in analytic_propagation_vec.m

[cfgMC] = initSimSSEM(cfgMC); %, Simulation, TLElaunchRepeat);   % Set population and launch rate


% PROPAGATOR
cfgMC.use_sgp4 = false;             % only 'false' is currently supported

% COLLISION
cfgMC.skipCollisions = 0;           % if 1, collision step is skipped in main_mc
cfgMC.max_frag = inf;

% Cube method Resolution level, note the statement below from J.C.Liou 2003
% "a dimension of 1% or less of the average semimajor axis of objects in the system is sufficient."
% 0.5% is used by "Assessing Collision Algorithms for the NewSpace Era"
%     cfgMC.CUBE_RES = mean(cfgMC.a_all(cfgMC.ap_all<(cfgMC.altitude_limit_up+radiusearthkm),1))*0.01;
cfgMC.CUBE_RES = 10;  % km
cfgMC.collision_alt_limit = 40000; %Ignoring satellites above 45000km for collision evaluation

% ATMOSPHERIC MODEL (pre-computed JB2008)
density_profile = 'static'; % The options are 'static' or 'JB2008'
cfgMC.density_profile = density_profile;
if strcmpi(cfgMC.density_profile,'JB2008')
    cfgMC = initJB2008(cfgMC);
end

% Save output file
cfgMC.save_output_file = 10;         % 0: don't save;
                                    % 1: save entire workspace
                                    % 2: sats_info, config, params
                                    % 3: summary (S_MC, N_MC, ...), config
                                    % 4: summary and collision stats (frag_info)
                                    % 5: just collision stats (frag_info)
                                    % 6: summary + average physical parameters
                                    % 7: Collision tracking (UROP)  
                                    % 10: entire workspace skipped every 10 timesteps - use main_mc2.m

cfgMC.saveMSnTimesteps = 73; % 73 = every 1 yr
                                    
paramSSEM.cd = 2.2;             % drag coefficient (flat plate model)
paramSSEM.rho0_re = 0.1570;     % B* constant (this number should already include the division by the earth radius)
cfgMC.paramSSEM = paramSSEM;    % update paramSSEM

tv = datestr(now, 'yyyymmddTHHMMSS');
filename_save = sprintf('Results/TLE_%s_%i.mat',tv,rngseed);
cfgMC.filename_save = filename_save;
cfgMC.n_save_checkpoint = 100; % save the results every n_save_checkpoint propagation step

return;

end


function cfgMCout = initJB2008(cfgMC)
% NB: it works just from March 2020 as a starting date of the propagation
load('dens_jb2008_032020_022224.mat')

dens_times=zeros(length(dens_highvar.month),1);

for k=1:length(dens_highvar.month)
    dens_times(k,1)= juliandate(datetime(dens_highvar.year(k),dens_highvar.month(k),0));
end

[dens_times2,alt2] = meshgrid(dens_times,dens_highvar.alt);

cfgMCout = cfgMC;
cfgMCout.param.dens_times = dens_times2;
cfgMCout.param.alt = alt2;
cfgMCout.param.dens_value = dens_highvar.dens;
end



function [cfgMCout,gm1,gm2,gm3] = initSimSSEM(cfg) % , Simulation, TLElaunchRepeat)
whichconst = 84;
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc( whichconst );

% TLEs-based initial population + random generation of the launch rate
% try
%     fn = which('initialized_01-2023.mat');  % 1 MB
%     % old version: initialized.mat; 25 MB
%     load(fn);  % just 'mat_sats'
%     fprintf('%i satellite entries loaded from %s\n', size(mat_sats,1),fn)
% 
%     % compute all Perigees and Apogees for filtering
%     a_all = mat_sats(:,1)*radiusearthkm;
%     e_all = mat_sats(:,2);
%     ap_all = a_all.*(1-e_all);
%     aa_all = a_all.*(1+e_all);
% catch
%     fprintf('current path: %s\n', pwd);
%     error('initial population (.mat) not found in path');
% end

getidx();

% Filter for desired altitude
% mat_sats(ap_all>(cfg.altitude_limit_up + radiusearthkm) | ap_all<(cfg.altitude_limit_low + radiusearthkm),:) = [];
% maxTLEssn = max(mat_sats(:,idx_ID));

% % remove RBs from iPop if flag is set
% if cfg.noRBs
%     mat_sats(mat_sats(:,idx_objectclass) == 5 | mat_sats(:,idx_objectclass)==6,:) = [];
%     if [size(cfg.paramSSEM.lam2,2) > 1 && sum(cfg.paramSSEM.lam2(:,2)) > 0]
%         error('lam2 has RB launches when cfg.noRBs flag is on!')
%     end
% end

% CREATE NEW OBJECTS (synthetic iPop)
% 1) From initialized.mat data, get physical param, and indexes by object type (g1,g2,g3)
% 2) Randomly select objects for OEs (correct num according to paramSSEM.lam2 and ipop)
% 3) Delete data to reset (SMA, radius, mass, obj number, controlled)
% 4) set data one by one

% [~,g1,g2,g3] = fillMassRadiusResample(mat_sats);   % Get indexes (payloads / RBs / debs) and GMs
% gm1 = g1.gm; gm2 = g2.gm; gm3 = g3.gm;

Smass = cfg.paramSSEM.mass(1);
Dmass = cfg.paramSSEM.mass(2);
Nmass = cfg.paramSSEM.mass(3);
Sradius = cfg.paramSSEM.radius(1);
Dradius = cfg.paramSSEM.radius(2);
Nradius = cfg.paramSSEM.radius(3);

nS = sum(cfg.paramSSEM.iPop(:,1));  % inputted initial population shape
nD = sum(cfg.paramSSEM.iPop(:,2));
nN = sum(cfg.paramSSEM.iPop(:,3));
% nB = sum(cfg.paramSSEM.iPop(:,4));

lS = sum(cfg.paramSSEM.lam2(:,1));  % inputted launch rates
% lB = sum(cfg.paramSSEM.lam2(:,2));

% matsatsSD = mat_sats(g1.allclass,:); 
% matsatsB = mat_sats(g2.allclass,:); 
% matsatsN = mat_sats(g3.allclass,:); 
% 
% selectS = matsatsSD(randi([1,size(matsatsSD,1)],nS,1),:); % select random sats from repeatLaunches
% selectD = matsatsSD(randi([1,size(matsatsSD,1)],nD,1),:);
% selectSlam = matsatsSD(randi([1,size(matsatsSD,1)],lS,1),:); % for lambda
% selectN = matsatsN(randi([1,size(matsatsN,1)],nN,1),:);

selectS = zeros(nS,24);
selectD = zeros(nD,24);
selectN = zeros(nN,24);
selectSlam = zeros(lS, 24);

% % Delete data to reset (SMA, radius, mass, obj number, controlled)
% selectS(:,[idx_a, idx_radius, idx_mass, idx_ID, idx_controlled]) = 0;
% selectD(:,[idx_a, idx_radius, idx_mass, idx_ID, idx_controlled]) = 0;
% selectN(:,[idx_a, idx_radius, idx_mass, idx_ID, idx_controlled]) = 0;
% selectSlam(:,[idx_a, idx_radius, idx_mass, idx_ID, idx_controlled]) = 0;

% scramble launch dates (random between epoch and epoch-missionLife)
% minDate =  jday(year(cfg.time0),month(cfg.time0),day(cfg.time0) - cfg.missionlifetime, 0,0,0);  % ~8 yrs prior to time0
minDate =  jday(year(cfg.time0) - cfg.missionlifetime, month(cfg.time0), day(cfg.time0), 0,0,0);  % ~8 yrs prior to time0; in JD
selectS(:,idx_launch_date) = minDate + rand(size(selectS,1),1) * cfg.YEAR2DAY * cfg.missionlifetime;
selectD(:,idx_launch_date) = minDate + rand(size(selectD,1),1) * cfg.YEAR2DAY * cfg.missionlifetime;
selectN(:,idx_launch_date) = minDate + rand(size(selectN,1),1) * cfg.YEAR2DAY * cfg.missionlifetime;
selectSlam(:,idx_launch_date) = minDate + rand(size(selectSlam,1),1) * cfg.YEAR2DAY;


% Select random SMA per shell
R02 = cfg.paramSSEM.R02;
re = cfg.paramSSEM.re;
defsC = ones(1,3);  % counter for iPop -- old: (1,4)
defsClam = ones(1,1);  % counter for lambda -- old: (1,2)

for ishell = 1 : numel(R02)-1
    defs = cfg.paramSSEM.iPop(ishell,:);
    defsLam = cfg.paramSSEM.lam2(ishell,:);
    if defs(1) > 0
        selectS(defsC(1):defsC(1) + defs(1) - 1, idx_a) = ...
            (rand(defs(1),1) * (R02(ishell+1) - R02(ishell)) + R02(ishell))./re + 1;
    end
    if defs(2) > 0
        selectD(defsC(2):defsC(2) + defs(2) - 1, idx_a) = ...
            (rand(defs(2),1) * (R02(ishell+1) - R02(ishell)) + R02(ishell))./re + 1;
    end
    if defs(3) > 0
        selectN(defsC(3):defsC(3) + defs(3) - 1, idx_a) = ...
            (rand(defs(3),1) * (R02(ishell+1) - R02(ishell)) + R02(ishell))./re + 1;
    end
    % if defs(4) > 0
    %     selectB(defsC(4):defsC(4) + defs(4) - 1, idx_a) = ...
    %         (rand(defs(4),1) * (R02(ishell+1) - R02(ishell)) + R02(ishell))./re + 1;
    % end
    if defsLam(1) > 0
        selectSlam(defsClam(1):defsClam(1) + defsLam(1) - 1, idx_a) = ...
            (rand(defsLam(1),1) * (R02(ishell+1) - R02(ishell)) + R02(ishell))./re + 1;
    end
    % if defsLam(2) > 0
    %     selectBlam(defsClam(2):defsClam(2) + defsLam(2) - 1, idx_a) = ...
    %         (rand(defsLam(2),1) * (R02(ishell+1) - R02(ishell)) + R02(ishell))./re + 1;
    % end
    defsC = defsC + defs;
    defsClam = defsClam + defsLam;
end

% Mark payloads as controlled
selectS(:,idx_controlled) = 1;
selectSlam(:,idx_controlled) = 1;
selectS(:,idx_a_desired) = selectS(:,idx_a);
selectSlam(:,idx_a_desired) = selectSlam(:,idx_a);

% mass
selectS(:,idx_mass) = Smass;
selectSlam(:,idx_mass) = Smass;
selectD(:,idx_mass) = Dmass;
selectN(:,idx_mass) = Nmass;

% radius
selectS(:,idx_radius) = Sradius;
selectSlam(:,idx_radius) = Sradius;
selectD(:,idx_radius) = Dradius;
selectN(:,idx_radius) = Nradius;

% missionlife
selectS(:,idx_missionlife) = cfg.missionlifetime;
selectSlam(:,idx_missionlife) = cfg.missionlifetime;

% obj type
selectS(:,idx_objectclass) = 1;
selectSlam(:,idx_objectclass) = 1;
selectD(:,idx_objectclass) = 1;
selectN(:,idx_objectclass) = 9;

% all together: ecco, inclo, bstar, mo, nodeo,  argpo, bstar, 
ipop = [selectS; selectD; selectN; selectSlam];  % include Lambda for now
ipop(:,idx_ecco) = 1e-6;
ipop(:,idx_inclo) = rand(size(ipop,1),1) * deg2rad(100);  % 0 to 100 deg
ipop(:,idx_mo) = rand(size(ipop,1),1) * 2*pi; 
ipop(:,idx_nodeo) = rand(size(ipop,1),1) * 2*pi; 
ipop(:,idx_argpo) = rand(size(ipop,1),1) * 2*pi; 
ipop(:,idx_bstar) = pi * 0.5 * 2.2 * ipop(:,idx_radius) .^2 ./ ipop(:,idx_mass) * 0.157;
ipop(:,idx_ID) = transpose(1:size(ipop,1));         % reset NORAD ID

lam2 = ipop(end-size(selectSlam,1)+1 : end,:);  % extract Slam
ipop(end-size(selectSlam,1)+1 : end,:) = []; 

    % 
    % constMatsat = nan(size(avglaunch,2), 24);
    % constMatsat(:,[idx_ecco, idx_inclo, idx_constel, idx_bstar, idx_mass, idx_radius, idx_missionlife]) = ...
    %     deg2rad(constellation.Inclination), constellation.constIndex, ...
    %     0.5 * 2.2 * constellation.radius .^2 ./ constellation.mass * 0.157, ...
    %     constellation.mass, constellation.radius, constellation.missionlife];
    % constMatsat(:, idx_a_desired) = constMatsat(:, idx_a);
    % constMatsat(:,[idx_controlled, idx_objectclass]) = [ ones(size(avglaunch))',  ones(size(avglaunch))'];
    % 
    % constLaunch = constMatsat(repelem(1:size(constMatsat,1), curConstLaunch), :  ); % populate constellation launch matsats
    % 
    % % bstar = 0.5 * 2.2 * ipop(:, idx_radius).^2 ./ ipop(:, idx_mass) * 0.157;
    % 
    % % Scramble launch time within this time step (dt_days), argpo, mo, nodeo
    % constLaunch(:,idx_launch_date) = jd + rand(size(constLaunch,1),1) * dt_days;
    % constLaunch(:,[idx_argpo, idx_mo, idx_nodeo]) = 2 * pi * rand(size(constLaunch,1),3);
    % 
    % param.maxID = param.maxID+size(out_future,1)+size(constLaunch,1);       %update maximum ID in population
    % count_tot_launches = count_tot_launches+size(out_future,1)+size(constLaunch,1);
    % out_future = [out_future; constLaunch];




% RESAMPLE physical param
% % ipop = multiplyInitialPop(cfg, [selectS; selectD; selectN; selectB], 1, g1,g2,g3);   % resampling method used, mult = x1
% % cfgMC.mat_sats = [selectS; selectD; selectN; selectB];
% ipop = fillMassRadiusResample([selectS; selectD; selectN; selectB],g1,g2,g3);
% 
% ipop(:, idx_a_desired) = ipop(:,idx_a);             % set current SMA as "desired SMA"
% ipop(:, idx_missionlife) = cfg.missionlifetime;     % set mission life


% RECALCULATE B* IF cfgMC.physicalBstar

% if cfg.physicalBstar
%     % Recalculate B* as 1/2 * Cd * A/m * rho0;  rho0 = 0.157 kg/m^2/RE
%     %   Physical def: C_0 = 1/2 * Cd * A/m * rho_0
%     %   TLE's use     C_0 = Bstar / 0.157 * rho;
%     % see line 89 in analytic_propagation_vec.m
%     bstari1 = ipop(:, idx_bstar);
%     bstarl1 = lam2(:, idx_bstar);
%     ipop(:, idx_bstar) = 0.5 * 2.2 * ipop(:, idx_radius).^2 ./ ipop(:, idx_mass) * 0.157;
%     lam2(:, idx_bstar) = 0.5 * 2.2 * lam2(:, idx_radius).^2 ./ lam2(:, idx_mass) * 0.157;    
%     bstari2 = ipop(:, idx_bstar);
%     bstarl2 = lam2(:, idx_bstar);
% end
%     figure; 
%     plot(bstari1,bstari2,'r-'); hold on; plot(bstarl1, bstarl2, 'b-'); legend('iPop','lambda');
%     xlabel('before'); ylabel('after');

%     figure;
%     subplot(211); histogram(bstari1,[-2:0.1:10],'edgealpha',0); a = gca; a.YScale = 'log'; title('B* Initial Pop');
%     hold on; histogram(bstari2,[-2:0.1:10],'edgealpha',0); yy = ylim; ylim([0.1,yy(2)]); legend('Before','After')
%     subplot(212); histogram(bstarl1,[-2:0.1:10],'edgealpha',0); hold on; histogram(bstarl2,[-2:0.1:10],'edgealpha',0); 
%     a = gca; a.YScale = 'log'; title('B* Launches'); yy = ylim; ylim([0.1,yy(2)]); legend('Before','After')


% figure;histogram(repLaunch(:,idx_objectclass));
% title(sprintf('Years %i to %i repeatedly launched', ...
%     cfg.launchRepeatYrs(1), cfg.launchRepeatYrs(2)));
% xlabel('Object type of launched objects');


% MAKE SURE IPOP AND LAM PERIGEES ARE WITHIN SPEC
% from prop_mit_vec: check_alt_ecc = out_mean_oe(:,1).*(1-out_mean_oe(:,2))>req+150 & out_mean_oe(:,2)<1; %check if decayed or hyperbolic
check_alt_ecc1 = find(ipop(:,idx_a).*(1-ipop(:,idx_ecco)) <= 1 + 150/radiusearthkm);
check_alt_ecc2 = find(ipop(:,2) >= 1); %check if decayed or hyperbolic
ipop(unique([check_alt_ecc1,check_alt_ecc2]),idx_ecco) = 0.001;  % circularize perigees < 150 km altitude, or hyperbolic objects

check_alt_ecc1 = find(lam2(:,idx_a).*(1-lam2(:,idx_ecco)) <= 1 + 150/radiusearthkm);
check_alt_ecc2 = find(lam2(:,2) >= 1);
lam2(unique([check_alt_ecc1,check_alt_ecc2]),idx_ecco) = 0.001;  % circularize


cfgMCout = cfg;
% cfgMCout.a_all = a_all;
% cfgMCout.ap_all = ap_all;
% cfgMCout.aa_all = aa_all;
cfgMCout.mat_sats = ipop;
cfgMCout.repeatLaunches = lam2;
% cfgMCout.additional_launches = additional_launches;
% cfgMCout.ind_launch = ind_launch;
% cfgMCout.ind_launch_add = ind_launch_add;
% cfgMCout.launch_model = launch_model;

end


function [outmatsat,g1,g2,g3] = fillMassRadiusResample(inmatsat,varargin)
    % Resample existing data to fill in missing info from DISCOS
    % Input: inmatsat -- matsat to be modified
    %        g1,g2,g3 -- optional; gX.gm defines the sampling dist (gmdistribution object)
    % object classes via objclass2int(1:12,2):
    %   P,PMRO,Pfrag,Pdeb,RB,RBMRO,RBfrag,RBdeb,Deb,OtherDeb,Unkwn,untracked
    
    % MATSATS DEFINITION
    getidx();
    
    % dobj = inmatsat(:,idx_objectclass);
    % for stats, combine into:
    %    payload {1}, RB {5}, debris {2,3,4,6,7,8,9,10,11}
    
    
    % fit data, save mean and covariances - if available, use given g1,g2,g3's GM model
    [g1,g2,g3] = getZeroGroups(inmatsat);
    if nargin > 1 && ~isempty(intersect(fieldnames(varargin{1}),'gm'))
        g1.gm = varargin{1}.gm;  % if supplied, use provided GM models
        g2.gm = varargin{2}.gm;
        g3.gm = varargin{3}.gm;
    else  % if no GM provided, calculate
        % g1.gm = gmdistribution(0,0);
        % g2.gm = gmdistribution(0,0);
        % g3.gm = gmdistribution(0,0);
        X = inmatsat(g1.nzno,[idx_radius, idx_mass]);
        if ~isempty(X)
            GMModel = fitgmdist(X,1);  g1.gm = GMModel;
        else
            g1.gm =  gmdistribution([],[]);
        end
        X = inmatsat(g2.nzno,[idx_radius, idx_mass]);
        if ~isempty(X)
            GMModel = fitgmdist(X,1);  g2.gm = GMModel;
        else
            g2.gm = gmdistribution([],[]);
        end
            X = inmatsat(g3.nzno,[idx_radius, idx_mass]);
        if ~isempty(X)
            GMModel = fitgmdist(X,1);  g3.gm = GMModel;
        else
            g3.gm = gmdistribution([],[]);
        end
    end
    
    % fill in empty data from mat_sats via sampling
    % -- current method: re-sample even if only one parameter is missing (radius or mass)
    % -- method: sample 2x N needed, remove any negative entries, randomly choose N
    outmatsat = inmatsat;
    if 0 ~= union(g1.zm, g1.zr)
        cursamp = random(g1.gm, numel(union(g1.zm, g1.zr))*2);
        cursamp(any(cursamp<=0,2),:) = []; % remove negative entries
        outmatsat(union(g1.zm, g1.zr),[idx_radius,idx_mass]) = cursamp(1:numel(union(g1.zm, g1.zr)),:);
    end
    if 0 ~= union(g2.zm, g2.zr)
        cursamp = random(g2.gm, numel(union(g2.zm, g2.zr))*2);
        cursamp(any(cursamp<=0,2),:) = []; % remove negative entries
        outmatsat(union(g2.zm, g2.zr),[idx_radius,idx_mass]) = cursamp(1:numel(union(g2.zm, g2.zr)),:);
    end
    if 0 ~= union(g3.zm, g3.zr)
        cursamp = random(g3.gm, numel(union(g3.zm, g3.zr))*2);
        cursamp(any(cursamp<=0,2),:) = []; % remove negative entries
        outmatsat(union(g3.zm, g3.zr),[idx_radius,idx_mass]) = cursamp(1:numel(union(g3.zm, g3.zr)),:);
    end
    fprintf('Resampled %i, %i, %i (g1,g2,g3) objects \n', numel(union(g1.zm, g1.zr)), ...
        numel(union(g2.zm, g2.zr)), numel(union(g3.zm, g3.zr)))
    
    return;

end

function [g1,g2,g3] = getZeroGroups(inmatsat)
getidx();

% mass vs radius
g1.allclass = []; % group 1: payloads; logical index
g2.allclass = []; % group 2: RBs
g3.allclass = []; % group 3: all debris

for ii = 1:12
    msinds = inmatsat(:,idx_objectclass) == ii; % all obj w/ objclass
    if ii == 1
        g1.allclass = find(msinds);                                 % all payload entries
        g1.zr = g1.allclass(inmatsat(g1.allclass,idx_radius) == 0); % index of zero radius
        g1.zm = g1.allclass(inmatsat(g1.allclass,idx_mass) == 0);   % index of zero mass
        g1.nz = setdiff(g1.allclass, union(g1.zr,g1.zm));           % index of non-zero radius and mass
        g1.nzno = g1.nz(~isoutlier(inmatsat(g1.nz,idx_radius)) ...  % non-zero, non-outlier
            & ~isoutlier(inmatsat(g1.nz,idx_mass)) );
    elseif ii == 5 | ii == 6
        g2.allclass = find(msinds);
        g2.zr = g2.allclass(inmatsat(g2.allclass,idx_radius) == 0); % all RB entries
        g2.zm = g2.allclass(inmatsat(g2.allclass,idx_mass) == 0);
        g2.nz = setdiff(g2.allclass, union(g2.zr,g2.zm));
        g2.nzno = g2.nz(~isoutlier(inmatsat(g2.nz,idx_radius)) ...
            & ~isoutlier(inmatsat(g2.nz,idx_mass)) );
    else
        g3.allclass = union(g3.allclass, find(msinds));             % all debris entries
    end
end
g3.zr = g3.allclass(inmatsat(g3.allclass,idx_radius) == 0);
g3.zm = g3.allclass(inmatsat(g3.allclass,idx_mass) == 0);
g3.nz = setdiff(g3.allclass, union(g3.zr, g3.zm));
g3.nzno = g3.nz(~isoutlier(inmatsat(g3.nz,idx_radius)) ...
    & ~isoutlier(inmatsat(g3.nz,idx_mass)) );
end

function mat_sats = multiplyInitialPop(cfgMC,inmatsat,mult,g1,g2,g3)
%     mult = cfgMC.xPop;
% randomly choose which satellites to add to initial population (per object type)
% note that physical parameters will be filled in by resampling (fillMassRadiusResample)
f1 = @(ac) randi(numel(ac), [floor(numel(ac) * mult),1]);     % repeats allowed
f2 = @(ac) randi(numel(ac), [floor(numel(ac) * (mult-1)),1]); % ensure all of original matsats are included
randinds = {};  acs = {[g1.allclass],[g2.allclass],[g3.allclass] };
for ind = 1:3
    if isempty(acs{ind})
        randinds{ind} = [];
    else
        if mult < 1
            randinds{ind} = f1(acs{ind});
        else
            randinds{ind} = f2(acs{ind});
        end
    end
end

if mult < 1  % simply select a subset of original matsats
    %         randind1 = randi(numel(g1.allclass), [floor(numel(g1.allclass) * mult),1]); % repeats allowed
    %         randind2 = randi(numel(g2.allclass), [floor(numel(g2.allclass) * mult),1]);
    %         randind3 = randi(numel(g3.allclass), [floor(numel(g3.allclass) * mult),1]);
    mat_sats = inmatsat([g1.allclass(randinds{1});g2.allclass(randinds{2});g3.allclass(randinds{3})],:);
else % if mult>1, ensure all of original matsats are included
    %         randind1 = randi(numel(g1.allclass), [floor(numel(g1.allclass) * (mult-1)),1]);
    %         randind2 = randi(numel(g2.allclass), [floor(numel(g2.allclass) * (mult-1)),1]);
    %         randind3 = randi(numel(g3.allclass), [floor(numel(g3.allclass) * (mult-1)),1]);
    % DEBUG
%     fprintf('max allclass g1/2/3: %i, %i, %i \n', ...
%         max(g1.allclass),max(g2.allclass),max(g3.allclass));
%     fprintf('size allclass of g1/2/3: %i, %i, %i \n', ...
%         size(g1.allclass,1),size(g2.allclass,1),size(g3.allclass,1));
%     fprintf('max randinds of 1/2/3: %i, %i, %i \n', ...
%         max(randinds{1}),max(randinds{2}),max(randinds{3}));
%     fprintf('max gX.allclass(randinds{X}) of 1/2/3: %i, %i, %i \n', ...
%         max(g1.allclass(randinds{1})), max(g2.allclass(randinds{2})), max(g3.allclass(randinds{3})));
%     fprintf(' Size of inmatsat: %i \n', ...
%         size(inmatsat,1));

    extrasats = inmatsat([g1.allclass(randinds{1});g2.allclass(randinds{2});g3.allclass(randinds{3})],:);
    % zero out mass and radius and fill in (resample)
    idx_mass = 8; idx_radius = 9;
    extrasats(:,[idx_mass, idx_radius]) = 0;
    if ~isempty(extrasats)
        if cfgMC.fillMassRadius > 0
            addmatsats = fillMassRadiusResample(extrasats,g1,g2,g3);
        else
            addmatsats = extrasats;
        end
        % scramble argperigee and mean motion
        idx_argpo = 5; idx_mo = 6;  % range: [0,2pi]
        addmatsats(:,[idx_argpo, idx_mo]) = 2 * pi * rand(size(addmatsats,1),2);
        mat_sats = [inmatsat; addmatsats];
    else
        mat_sats = inmatsat;  % skip addmatsats if extrasats is empty
    end
end
end

function lam = getLamFromMatsats(cfg)
% cfg fields used: launchRepeatYrs, R02, repeatLaunches, paramSSEM
inmatsat = cfg.repeatLaunches;
getidx();

if cfg.launchRepeatSmooth  % just one year's worth of launches in cfg.repeatLaunches
    lam = transpose(histcounts(inmatsat(inmatsat(:,idx_controlled) == 1, idx_a), ...
        cfg.paramSSEM.R02 / cfg.paramSSEM.re + 1));
    % output a column vector
else
    lam = transpose(histcounts(inmatsat(inmatsat(:,idx_controlled) == 1, idx_a), ...
        cfg.paramSSEM.R02 / cfg.paramSSEM.re + 1) / (range(cfg.launchRepeatYrs) + 1));
end

end
