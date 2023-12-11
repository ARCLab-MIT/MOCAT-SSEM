function asem = matsats2ASEM(matsats, idxs)
% matsats to ASEM matrix

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

getidx();
radiusearthkm = 6378.137;

% idx_a = 1; idx_ecco = 2; idx_inclo = 3; idx_nodeo = 4; idx_argpo = 5; idx_mo = 6; idx_bstar = 7; idx_mass = 8; idx_radius = 9;
%     idx_error = 10; idx_controlled = 11; idx_a_desired = 12; idx_missionlife = 13; idx_constel = 14; idx_date_created = 15; idx_launch_date = 16;
%     idx_r = [17 18 19]; idx_v = [20 21 22]; idx_objectclass = 23; idx_ID = 24;
asem = [];
asem(:,3) = matsats(:,idx_a) .* radiusearthkm;  % SMA [km]
asem(:,4) = matsats(:,idx_ecco);                % Ecc [-]
asem(:,5) = rad2deg(matsats(:,idx_inclo));      % Inc [rad]
asem(:,6) = rad2deg(matsats(:,idx_nodeo));      % RAAN [rad]
asem(:,7) = rad2deg(matsats(:,idx_argpo));      % AOP [rad]
asem(:,8) = rad2deg(matsats(:,idx_mo));         % Mean Anomaly [rad]
asem(:,14) = matsats(:,idx_mass);               % Mass [kg]
asem(:,15) = matsats(:,idx_radius)*2;           % Diameter [m]
%asem(:,13) = matsats(:,idx_bstar) .* asem(:,14) / (1/2 * 2.2 *0.157e6); % area [m^2] derived from bstar
%     % calculate B* as Bstar = 1/2 * Cd * A/m * 0.157e6
asem(:,13) = pi*matsats(:,idx_radius).^2; % area [m^2] derived from bstar

asem(matsats(:, idx_objectclass) == 5, 10) = 1; % RB; see objclass2int
asem(matsats(:, idx_objectclass) == 1, 10) = 2; % payload

deb_spec = [2:4,6:12];
debidx = ismember(matsats(:,idx_objectclass), deb_spec);
asem(debidx,10) = 3;                            % Deb


asem(:,12) = 0;
asem(matsats(:,idx_controlled) == 1, 12) = 3;   % Controlled
asem(:,2) = matsats(:,idx_launch_date);         % start or launch date [JD]
asem(:,1) = matsats(:,idx_ID); 
asem(:,9) = asem(:,2) + matsats(:,idx_missionlife)*365.2425;  % [JD]

asem(:,16)  = 1;

cols = cellstr(["adept_id", "epoch_start", "sma", "ecc", "inc", "raan", "aop", "ma", "epoch_end", "obj_type", "disp_option", "stationkeeping", "area", "mass", "size", "weight"]);
asem = array2table(asem, 'VariableNames', cols);

end