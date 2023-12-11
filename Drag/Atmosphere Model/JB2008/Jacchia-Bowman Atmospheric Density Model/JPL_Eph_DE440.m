%--------------------------------------------------------------------------
%
% JPL_Eph_DE436: Computes the sun, moon, and nine major planets' equatorial
%                position using JPL Ephemerides
%
% Inputs:
%   Mjd_TDB         Modified julian date of TDB
%
% Output:
%   r_Earth,r_SunSSB(solar system barycenter (SSB)),r_Mercury,r_Venus,
%   r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,
%   r_Sun(geocentric equatorial position ([m]) referred to the
%   International Celestial Reference Frame (ICRF))
%
% Last modified:   2022/05/18   Meysam Mahooti
% 
%--------------------------------------------------------------------------
function [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
          r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE440(Mjd_TDB)

global PC

JD = Mjd_TDB + 2400000.5;

i = find(PC(:,1)<=JD & JD<=PC(:,2),1,'first'); 
PCtemp = PC(i,:);

t1 = PCtemp(1)-2400000.5; % MJD at start of interval

dt = Mjd_TDB - t1;

temp = (231:13:270);
Cx_EarthMoon = PCtemp(temp(1):temp(2)-1);
Cy_EarthMoon = PCtemp(temp(2):temp(3)-1);
Cz_EarthMoon = PCtemp(temp(3):temp(4)-1);
temp = temp+39;
Cx = PCtemp(temp(1):temp(2)-1);
Cy = PCtemp(temp(2):temp(3)-1);
Cz = PCtemp(temp(3):temp(4)-1);
Cx_EarthMoon = [Cx_EarthMoon,Cx];
Cy_EarthMoon = [Cy_EarthMoon,Cy];
Cz_EarthMoon = [Cz_EarthMoon,Cz];    
if (0<=dt && dt<=16)
    j=0;
    Mjd0 = t1;
elseif(16<dt && dt<=32)
    j=1;
    Mjd0 = t1+16*j;
end
r_EarthMoon = 1e3*Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, Cx_EarthMoon(13*j+1:13*j+13),...
                     Cy_EarthMoon(13*j+1:13*j+13), Cz_EarthMoon(13*j+1:13*j+13))';

temp = (441:13:480);
Cx_Moon = PCtemp(temp(1):temp(2)-1);
Cy_Moon = PCtemp(temp(2):temp(3)-1);
Cz_Moon = PCtemp(temp(3):temp(4)-1);
for i=1:7
    temp = temp+39;
    Cx = PCtemp(temp(1):temp(2)-1);
    Cy = PCtemp(temp(2):temp(3)-1);
    Cz = PCtemp(temp(3):temp(4)-1);   
    Cx_Moon = [Cx_Moon,Cx];
    Cy_Moon = [Cy_Moon,Cy];
    Cz_Moon = [Cz_Moon,Cz];    
end
if (0<=dt && dt<=4)
    j=0;
    Mjd0 = t1;
elseif(4<dt && dt<=8)
    j=1;
    Mjd0 = t1+4*j;
elseif(8<dt && dt<=12)
    j=2;
    Mjd0 = t1+4*j;
elseif(12<dt && dt<=16)
    j=3;
    Mjd0 = t1+4*j;
elseif(16<dt && dt<=20)
    j=4;
    Mjd0 = t1+4*j;
elseif(20<dt && dt<=24)
    j=5;
    Mjd0 = t1+4*j;
elseif(24<dt && dt<=28)
    j=6;
    Mjd0 = t1+4*j;
elseif(28<dt && dt<=32)
    j=7;
    Mjd0 = t1+4*j;
end
r_Moon = 1e3*Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, Cx_Moon(13*j+1:13*j+13),...
                    Cy_Moon(13*j+1:13*j+13), Cz_Moon(13*j+1:13*j+13))';

temp = (753:11:786);
Cx_Sun = PCtemp(temp(1):temp(2)-1);
Cy_Sun = PCtemp(temp(2):temp(3)-1);
Cz_Sun = PCtemp(temp(3):temp(4)-1);
temp = temp+33;
Cx = PCtemp(temp(1):temp(2)-1);
Cy = PCtemp(temp(2):temp(3)-1);
Cz = PCtemp(temp(3):temp(4)-1);   
Cx_Sun = [Cx_Sun,Cx];
Cy_Sun = [Cy_Sun,Cy];
Cz_Sun = [Cz_Sun,Cz];
if (0<=dt && dt<=16)
    j=0;
    Mjd0 = t1;
elseif(16<dt && dt<=32)
    j=1;
    Mjd0 = t1+16*j;
end
r_Sun = 1e3*Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, Cx_Sun(11*j+1:11*j+11),...
                   Cy_Sun(11*j+1:11*j+11), Cz_Sun(11*j+1:11*j+11))';

temp = (3:14:45);
Cx_Mercury = PCtemp(temp(1):temp(2)-1);
Cy_Mercury = PCtemp(temp(2):temp(3)-1);
Cz_Mercury = PCtemp(temp(3):temp(4)-1);
for i=1:3
    temp = temp+42;
    Cx = PCtemp(temp(1):temp(2)-1);
    Cy = PCtemp(temp(2):temp(3)-1);
    Cz = PCtemp(temp(3):temp(4)-1);
    Cx_Mercury = [Cx_Mercury,Cx];
    Cy_Mercury = [Cy_Mercury,Cy];
    Cz_Mercury = [Cz_Mercury,Cz];    
end
if (0<=dt && dt<=8)
    j=0;
    Mjd0 = t1;
elseif(8<dt && dt<=16)
    j=1;
    Mjd0 = t1+8*j;
elseif (16<dt && dt<=24)
    j=2;
    Mjd0 = t1+8*j;
elseif(24<dt && dt<=32)
    j=3;
    Mjd0 = t1+8*j;
end
r_Mercury = 1e3*Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, Cx_Mercury(14*j+1:14*j+14),...
                       Cy_Mercury(14*j+1:14*j+14), Cz_Mercury(14*j+1:14*j+14))';

temp = (171:10:201);
Cx_Venus = PCtemp(temp(1):temp(2)-1);
Cy_Venus = PCtemp(temp(2):temp(3)-1);
Cz_Venus = PCtemp(temp(3):temp(4)-1);
temp = temp+30;
Cx = PCtemp(temp(1):temp(2)-1);
Cy = PCtemp(temp(2):temp(3)-1);
Cz = PCtemp(temp(3):temp(4)-1);
Cx_Venus = [Cx_Venus,Cx];
Cy_Venus = [Cy_Venus,Cy];
Cz_Venus = [Cz_Venus,Cz];
if (0<=dt && dt<=16)
    j=0;
    Mjd0 = t1;
elseif(16<dt && dt<=32)
    j=1;
    Mjd0 = t1+16*j;
end
r_Venus = 1e3*Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, Cx_Venus(10*j+1:10*j+10),...
                     Cy_Venus(10*j+1:10*j+10), Cz_Venus(10*j+1:10*j+10))';

temp = (309:11:342);
Cx_Mars = PCtemp(temp(1):temp(2)-1);
Cy_Mars = PCtemp(temp(2):temp(3)-1);
Cz_Mars = PCtemp(temp(3):temp(4)-1);
j=0;
Mjd0 = t1;
r_Mars = 1e3*Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, Cx_Mars(11*j+1:11*j+11),...
                    Cy_Mars(11*j+1:11*j+11), Cz_Mars(11*j+1:11*j+11))';

temp = (342:8:366);
Cx_Jupiter = PCtemp(temp(1):temp(2)-1);
Cy_Jupiter = PCtemp(temp(2):temp(3)-1);
Cz_Jupiter = PCtemp(temp(3):temp(4)-1);
j=0;
Mjd0 = t1;
r_Jupiter = 1e3*Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32, Cx_Jupiter(8*j+1:8*j+8),...
                       Cy_Jupiter(8*j+1:8*j+8), Cz_Jupiter(8*j+1:8*j+8))';

temp = (366:7:387);
Cx_Saturn = PCtemp(temp(1):temp(2)-1);
Cy_Saturn = PCtemp(temp(2):temp(3)-1);
Cz_Saturn = PCtemp(temp(3):temp(4)-1);
j=0;
Mjd0 = t1;
r_Saturn = 1e3*Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, Cx_Saturn(7*j+1:7*j+7),...
                      Cy_Saturn(7*j+1:7*j+7), Cz_Saturn(7*j+1:7*j+7))';

temp = (387:6:405);
Cx_Uranus = PCtemp(temp(1):temp(2)-1);
Cy_Uranus = PCtemp(temp(2):temp(3)-1);
Cz_Uranus = PCtemp(temp(3):temp(4)-1);
j=0;
Mjd0 = t1;
r_Uranus = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Uranus(6*j+1:6*j+6),...
                      Cy_Uranus(6*j+1:6*j+6), Cz_Uranus(6*j+1:6*j+6))';

temp = (405:6:423);
Cx_Neptune = PCtemp(temp(1):temp(2)-1);
Cy_Neptune = PCtemp(temp(2):temp(3)-1);
Cz_Neptune = PCtemp(temp(3):temp(4)-1);
j=0;
Mjd0 = t1;
r_Neptune = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Neptune(6*j+1:6*j+6),...
                       Cy_Neptune(6*j+1:6*j+6), Cz_Neptune(6*j+1:6*j+6))';

temp = (423:6:441);
Cx_Pluto = PCtemp(temp(1):temp(2)-1);
Cy_Pluto = PCtemp(temp(2):temp(3)-1);
Cz_Pluto = PCtemp(temp(3):temp(4)-1);
j=0;
Mjd0 = t1;
r_Pluto = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Pluto(6*j+1:6*j+6),...
                     Cy_Pluto(6*j+1:6*j+6), Cz_Pluto(6*j+1:6*j+6))';

temp = (819:10:839);
Cx_Nutations = PCtemp(temp(1):temp(2)-1);
Cy_Nutations = PCtemp(temp(2):temp(3)-1);
for i=1:3
    temp = temp+20;
    Cx = PCtemp(temp(1):temp(2)-1);
    Cy = PCtemp(temp(2):temp(3)-1);
    Cx_Nutations = [Cx_Nutations,Cx];
    Cy_Nutations = [Cy_Nutations,Cy];
end
if (0<=dt && dt<=8)
    j=0;
    Mjd0 = t1;
elseif(8<dt && dt<=16)
    j=1;
    Mjd0 = t1+8*j;
elseif (16<dt && dt<=24)
    j=2;
    Mjd0 = t1+8*j;
elseif(24<dt && dt<=32)
    j=3;
    Mjd0 = t1+8*j;
end
Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Nutations(10*j+1:10*j+10),...
                   Cy_Nutations(10*j+1:10*j+10),zeros(10,1))';

temp = (899:10:929);
Cx_Librations = PCtemp(temp(1):temp(2)-1);
Cy_Librations = PCtemp(temp(2):temp(3)-1);
Cz_Librations = PCtemp(temp(3):temp(4)-1);
for i=1:3
    temp = temp+30;
    Cx = PCtemp(temp(1):temp(2)-1);
    Cy = PCtemp(temp(2):temp(3)-1);
    Cz = PCtemp(temp(3):temp(4)-1);
    Cx_Librations = [Cx_Librations,Cx];
    Cy_Librations = [Cy_Librations,Cy];
    Cz_Librations = [Cz_Librations,Cz];    
end
if (0<=dt && dt<=8)
    j=0;
    Mjd0 = t1;
elseif(8<dt && dt<=16)
    j=1;
    Mjd0 = t1+8*j;
elseif (16<dt && dt<=24)
    j=2;
    Mjd0 = t1+8*j;
elseif(24<dt && dt<=32)
    j=3;
    Mjd0 = t1+8*j;
end
Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Librations(10*j+1:10*j+10),...
                    Cy_Librations(10*j+1:10*j+10), Cz_Librations(10*j+1:10*j+10))';
EMRAT = 81.3005682214972154; % DE440
EMRAT1 = 1/(1+EMRAT);
r_Earth = r_EarthMoon-EMRAT1*r_Moon;
r_Mercury = -r_Earth+r_Mercury;
r_Venus = -r_Earth+r_Venus;
r_Mars = -r_Earth+r_Mars;
r_Jupiter = -r_Earth+r_Jupiter;
r_Saturn = -r_Earth+r_Saturn;
r_Uranus = -r_Earth+r_Uranus;
r_Neptune = -r_Earth+r_Neptune;
r_Pluto = -r_Earth+r_Pluto;
r_SunSSB = r_Sun;
r_Sun = -r_Earth+r_Sun;

end

