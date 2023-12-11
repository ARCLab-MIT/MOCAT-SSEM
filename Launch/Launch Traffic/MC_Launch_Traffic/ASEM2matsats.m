function matsats = ASEM2matsats(fname)
% fname = 'start_1cm_to_10cm.asem'


% ASEM Start file format has 16 columns
%   COL#    DESCRIPTON
%   ----    ----------
%   [1]     Object ID
%   [2]     Julian Date (JD) Start Epoch
%   [3:8]   SMA, ECC, INC, RAAN, AOP, MA [km, deg]
%   [9]     Julian Date (JD) Final Epoch
%   [10]     Object type flag (1=roc,2=sat,3=deb,4=??)
%   [11]    Disposal option flag (0 for popZERO)
%   [12]    Stationkeeping flag
%             (0 = none, 1 = long only, 2 = long & inc, 3 = turn off drag)
%   [13]    Area (m^2)
%   [14]    Mass (kg)
%   [15]    Diameter (m)
%   [16]    Weighting Factor

try
    asem = readmatrix(fname,'FileType','text');
catch
    error('Cannot import %s via readmatrix',fname);
end

getidx();
radiusearthkm = 6378.137;

% idx_a = 1; idx_ecco = 2; idx_inclo = 3; idx_nodeo = 4; idx_argpo = 5; idx_mo = 6; idx_bstar = 7; idx_mass = 8; idx_radius = 9;
%     idx_error = 10; idx_controlled = 11; idx_a_desired = 12; idx_missionlife = 13; idx_constel = 14; idx_date_created = 15; idx_launch_date = 16;
%     idx_r = [17 18 19]; idx_v = [20 21 22]; idx_objectclass = 23; idx_ID = 24;
matsats = [];
matsats(:,idx_a) = asem(:,3)./radiusearthkm;    % SMA [Earth radii]
matsats(:,idx_ecco) = asem(:,4);                % Ecc [-]
matsats(:,idx_inclo) = deg2rad(asem(:,5));      % Inc [rad]
matsats(:,idx_nodeo) = deg2rad(asem(:,6));       % RAAN [rad]
matsats(:,idx_argpo) = deg2rad(asem(:,7));       % AOP [rad]
matsats(:,idx_mo) = deg2rad(asem(:,8));          % Mean anomaly [rad]
matsats(:,idx_bstar) = 1/2 * 2.2 * asem(:,13)./asem(:,14)*0.157e6;
                                                % Bstar [1/earthradii]
     % calculate B* as Bstar = 1/2 * Cd * A/m * 0.157e6
matsats(:,idx_mass) = asem(:,14);               % Mass [kg]
matsats(:,idx_radius) = asem(:,15)/2;           % Radius [m]
matsats(:,idx_error) = 0;

% Below: see objclass2int
matsats(asem(:,10) == 1, idx_objectclass) = 5;  % RB
matsats(asem(:,10) == 2, idx_objectclass) = 1;  % Payloads
matsats(asem(:,10) == 3, idx_objectclass) = 10; % Debris
matsats(asem(:,10) == 4, idx_objectclass) = 11; % unknown

matsats(:,idx_controlled) = ~asem(:,12);
matsats(:,idx_a_desired) = matsats(:,idx_a);  
matsats(:,idx_missionlife) = (asem(:,9) - asem(:,2))/365.2425;  % [yrs]
matsats(:,idx_constel) = 0;                     % constellation flag
matsats(:,idx_date_created) = asem(:,2);        % [JD]
matsats(:,idx_launch_date) = asem(:,2);         % [JD]
matsats(:,idx_r) = 0;
matsats(:,idx_v) = 0;
matsats(:,idx_ID) = asem(:,1);

end