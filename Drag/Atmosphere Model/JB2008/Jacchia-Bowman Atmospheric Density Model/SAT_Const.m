%--------------------------------------------------------------------------
%
% SAT_Const: Definition of astronomical and mathematical constants
%
% Last modified:   2022/05/18   Meysam Mahooti
%
%--------------------------------------------------------------------------

% Mathematical constants
const.pi2       = 2*pi;                % 2pi
const.Rad       = pi/180;              % Radians per degree
const.Deg       = 180/pi;              % Degrees per radian
const.Arcs      = 3600*180/pi;         % Arcseconds per radian

% General
const.MJD_J2000 = 51544.5;             % Modified Julian Date of J2000
const.T_B1950   = -0.500002108;        % Epoch B1950
const.c_light   = 299792457.999999984; % Speed of light  [m/s]; DE440
const.AU        = 149597870699.999988; % Astronomical unit [m]; DE440

% Physical parameters of the Earth, Sun and Moon

% Equatorial radius and flattening
const.R_Earth   = 6378.137e3;      % Earth's radius [m]; WGS-84
const.f_Earth   = 1/298.257223563; % Flattening; WGS-84
const.R_Sun     = 696000.0e3;      % Sun's radius [m]; DE440
const.R_Moon    = 1738.0e3;        % Moon's radius [m]; DE440

% Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
const.omega_Earth = 15.04106717866910/3600*const.Rad; % [rad/s]; WGS-84

% Gravitational coefficients
const.GM_Earth   = 398600.4415e9;         			   % [m^3/s^2]; GGM03C & GGM03S
const.GM_Sun     = 132712440041.279419e9; 			   % [m^3/s^2]; DE440
const.GM_Moon    = const.GM_Earth/81.3005682214972154; % [m^3/s^2]; DE440
const.GM_Mercury = 22031.868551e9; 		  			   % [m^3/s^2]; DE440
const.GM_Venus   = 324858.592000e9;       			   % [m^3/s^2]; DE440
const.GM_Mars    = 42828.375816e9;	      			   % [m^3/s^2]; DE440
const.GM_Jupiter = 126712764.100000e9;    			   % [m^3/s^2]; DE440
const.GM_Saturn  = 37940584.841800e9;     			   % [m^3/s^2]; DE440
const.GM_Uranus  = 5794556.400000e9;      			   % [m^3/s^2]; DE440
const.GM_Neptune = 6836527.100580e9;      			   % [m^3/s^2]; DE440
const.GM_Pluto   = 975.500000e9;	      			   % [m^3/s^2]; DE440

% Solar radiation pressure at 1 AU
const.P_Sol = 1367/const.c_light; % [N/m^2] (1367 W/m^2); IERS

