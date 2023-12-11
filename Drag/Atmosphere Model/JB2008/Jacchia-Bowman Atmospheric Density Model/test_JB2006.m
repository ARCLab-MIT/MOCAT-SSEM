%--------------------------------------------------------------------------
%
%     TEST CASE FOR RUNNING JB2006
% 
% Last modified:   2022/09/11   Meysam Mahooti
%
%--------------------------------------------------------------------------
clc
clear
format long g

global const PC

SAT_Const
constants
load DE440Coeff.mat
PC = DE440Coeff;

% read Earth orientation parameters
fid = fopen('eop19620101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

% read solar storm indices
fid = fopen('SOLFSMY.txt','r');
%  ---------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c Ssrc
%  ---------------------------------------------------------------------------
SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);

% read Ap data
fid = fopen('SOLRESAP.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD  F10 F10B Ap1 to Ap8
%  ------------------------------------------------------------------------
APdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[12 inf]);
fclose(fid);

year = 2001;
doy = 200;
[month,day,hour,minute,sec] = days2mdh(year,doy);
MJD = Mjday(year,month,day,hour,minute,sec);

% SET SOLAR INDICES
JD = floor(MJD+2400000.5);
i = find(JD==SOLdata(3,:),1,'first');
SOL = SOLdata(:,i);
F10B = SOL(5);
S10B = SOL(7);
XM10 = SOL(8);

% USE 1 DAY LAG FOR EUV
JD = floor(MJD-1+2400000.5);
i = find(JD==SOLdata(3,:),1,'first');
SOL = SOLdata(:,i);
F10 = SOL(4);
S10 = SOL(6);

% USE 5 DAY LAG FOR MG FUV INFLUENCE
SOL = SOLdata(:,i-4);
XM10B = SOL(9);

% USE 6.7 HR LAG FOR Ap INFLUENCE
[year,month,day,hour,minute,sec] = invjday(MJD-6.7/24+2400000.5);
doy = finddays(year,month,day,hour,minute,sec);
i = find(year==APdata(1,:) & floor(doy)==APdata(2,:),1,'first');
AP = APdata(5:12,i);

if hour < 3
    Ap = AP(1);
elseif hour < 6
    Ap = AP(2);
elseif hour < 9
    Ap = AP(3);
elseif hour < 12
    Ap = AP(4);
elseif hour < 15
    Ap = AP(5);
elseif hour < 18
    Ap = AP(6);
elseif hour < 21
    Ap = AP(7);
else
    Ap = AP(8);
end

GEO(1) = F10;
GEO(2) = F10B;
GEO(3) = Ap;

% SET POINT OF INTEREST LOCATION (RADIANS AND KM)
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
year = 2001;
doy = 200;
[month,day,hour,minute,sec] = days2mdh(year,doy);
[DJMJD0, DATE] = iauCal2jd(year, month, day);
TIME = (60*(60*hour+minute)+sec)/86400;
UTC = DATE+TIME;
TT = UTC+TT_UTC/86400;
TUT = TIME+UT1_UTC/86400;
UT1 = DATE+TUT;
GWRAS = iauGmst06(DJMJD0, UT1, DJMJD0, TT);
XLON = 60*const.Rad;
SAT(1) = mod(GWRAS + XLON, 2*pi);
SAT(2) = -70*const.Rad; 
SAT(3) = 400;

% Difference between ephemeris time and universal time
% JD = MJD+2400000.5;
% [year, month, day, hour, minute, sec] = invjday(JD);
% days = finddays(year, month, day, hour, minute, sec);
% ET_UT = ETminUT(year+days/365.25);
% MJD_ET = MJD+ET_UT/86400;
% [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
%  r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE440(MJD_ET);

MJD_TDB = Mjday_TDB(TT);
[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
 r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE440(MJD_TDB);

% SET Sun's right ascension and declination (RADIANS)
ra_Sun  = atan2(r_Sun(2), r_Sun(1));
dec_Sun = atan2(r_Sun(3), sqrt(r_Sun(1)^2+r_Sun(2)^2));
SUN(1)  = ra_Sun;
SUN(2)  = dec_Sun;

% COMPUTE DENSITY KG/M3 RHO
[TEMP,RHO] = JB2006(MJD,SUN,SAT,GEO,S10,S10B,XM10,XM10B)

