%***********************************************************************
%     Jacchia-Bowman 2006 Atmospheric Model
%
%     This is the CIRA "Integration Form" of a Jacchia Model.
%     There are no tabular values of density.  Instead, the barometric
%     equation and diffusion equation are integrated numerically using
%     the Newton-Coates method to produce the density profile up to the
%     input position.
%
%     INPUT:
%
%           MJD    : Date and Time, in modified Julian Days
%                    and Fraction (MJD = JD-2400000.5)
%           SUN(1) : Right Ascension of Sun (radians)
%           SUN(2) : Declination of Sun (radians)
%           SAT(1) : Right Ascension of Position (radians)
%           SAT(2) : Geocentric Latitude of Position (radians)
%           SAT(3) : Height of Position (km)
%           GEO(1) : 10.7-cm Solar Flux (1.0E-22*Watt/(M^2*Hertz))
%                    (Tabular time 1.0 day earlier)
%           GEO(2) : 10.7-cm Solar Flux, ave.
%                    81-day centered on the input time
%           GEO(3) : Geomagnetic Planetary 3-hour Index
%                    A-sub-p for a tabular time 0.279 days earlier
%                    (6.7 hours earlier)
%           S10    : EUV index (26-34 nm) scaled to F10
%                    (Tabular time 1.0 day earlier)
%           S10B   : EUV 81-day ave. centered index
%           XM10   : MG2 index scaled to F10
%           XM10B  : MG2 81-day ave. centered index
%                    (Tabular time 5.0 days earlier)
%
%     OUTPUT:
%
%           TEMP(1): Exospheric Temperature above Input Position (deg K)
%           TEMP(2): Temperature at Input Position (deg K)
%           RHO    : Total Mass-Desnity at Input Position (kg/m^3)
%
%
%     JB2006 Model Development: (Ref. 1)
%
%      1. Start with the CIRA72 model (Jacchia 71).
%
%      2. Convert to a CIRA70 model replacing equations from Jacchia 70
%         model (Ref. 5)
%
%      3. Replace Tc equation using new solar indices (Ref. 1 and 2)
%
%      4. Replace semiannual equation with new global model based
%         on F10B (Ref. 1 and 3)
%
%      5. Add correction for local solar time and latitude errors
%         (Ref. 1)
%         Added smooth transition between altitude bands
%
%      6. Add high altitude ( z > 1500 km ) correction
%         (Ref. 1 and 4)
%
%
%         All equation references below refer to the original
%         Jacchia model papers.
%
%
%     References:
%
%      1. Bowman, Bruce R., etc. : "A New Empirical Thermospheric
%         Density Model JB2006 Using New Solar Indices",
%         AIAA/AAS Astrodynamics Specialists Conference, Keystone, CO,
%         21-24 Aug 2006, (Paper AIAA 2006-6166).
%
%      2. Bowman, Bruce R., etc. : "Improvements in Modeling
%         Thermospheric Densities Using New EUV and FUV Solar Indices",
%         AAS/AIAA Space Flight Mechanics Meeting, Tampa, FL,
%         23-26 Jan 2006, (Paper AAS 06-237).
%
%      3. Bowman, Bruce R.: "The Semiannual Thermospheric Density
%         Variation From 1970 to 2002 Between 200-1100 km",
%         AAS/AIAA Space Flight Mechanics Meeting, Maui, HI,
%         8-12 Feb 2004, (Paper AAS 04-174).
%
%      4. Bowman, Bruce R.; "Atmospheric Density Variations at
%         1500 km to 4000 km Height Determined from Long Term
%         Orbit Perturbation Analysis", AAS/AIAA Space Flight
%         Mechanics Meeting, Santa Barbara, CA, 11-14 Feb 2001,
%         (Paper AAS 01-132).
%
%      5. Jacchia, Luigi G.; "New Static Models of the
%         Thermosphere and Exosphere with Empirical Temperature
%         Profiles", (Smithsonian Astrophysical Observatory
%         Special Report 313), 6 May 1970.
%
%
%     Written by: Bruce R Bowman (HQ AFSPC, Space Analysis Division),
%                 Feb 2006
%
%         REV A - Oct 2006
%                 Smoothing of density corrections and scale height
%                 through different altitude bands in the latitude-
%                 local time correction subroutine DTSUB
%                 dTx correction replaced with dTc correction
%
% Last modified:   2018/01/27   M. Mahooti
%
%***********************************************************************
function [TEMP,RHO] = JB2006 (MJD,SUN,SAT,GEO,S10,S10B,XM10,XM10B)

% The alpha are the thermal diffusion coefficients in Eq. (6)
ALPHA = [0,0,0,0,-0.38];

% AL10 is log(10.0)
AL10 = 2.3025851;

% The AMW are the molecular weights in order: N2, O2, O, Ar, He & H
AMW = [28.0134,31.9988,15.9994,39.9480,4.0026,1.00797];

% AVOGAD is Avogadro's number in mks units (molecules/kmol)
AVOGAD = 6.02257e26;

TWOPI  = 2*pi;
PIOV2  = 1.5707963;

% The FRAC are the assumed sea-level volume fractions in order:
% N2, O2, Ar, and He
FRAC = [0.78110,0.20955,9.3400e-3,1.2890e-5];

% RSTAR is the universal gas-constant in mks units (joules/K/kmol)
RSTAR = 8314.32;

% The R# are values used to establish height step sizes in
% the regimes 90km to 105km, 105km to 500km and 500km upward.
R1 = 0.010;
R2 = 0.025;
R3 = 0.075;

% The WT are weights for the Newton-Cotes Five-Point Quad. formula
WT =[0.311111111111111,1.422222222222222,0.533333333333333,...
     1.422222222222222,0.311111111111111];

% The CHT are coefficients for high altitude density correction
CHT = [0.22,-0.20e-02,0.115e-02,-0.211e-05];
DEGRAD = pi/180;

F10B = GEO(2);
AP   = GEO(3);

% Equation (14)
TSUBC = 379. + 3.353*GEO(2) + 0.358*(GEO(1)-GEO(2)) + ...
        2.094*(S10-S10B) + 0.343*(XM10-XM10B);

% Equation (15)
ETA =   0.5 * abs(SAT(2) - SUN(2));
THETA = 0.5 * abs(SAT(2) + SUN(2));

% Equation (16)
H = SAT(1) - SUN(1);
TAU = H - 0.64577182 + 0.10471976 * sin(H + 0.75049158);
GLAT  = SAT(2);
ZHT   = SAT(3);
GLST  = H + pi;
GLSTHR = (GLST/DEGRAD)*(24/360);

if (GLSTHR>=24)
    GLSTHR = GLSTHR - 24;
end
if (GLSTHR< 0)
    GLSTHR = GLSTHR + 24;
end

% Equation (17)
C = cos(ETA)^2.5;
S = sin(THETA)^2.5;

DF = S + (C - S) * abs(cos(0.5 * TAU))^3;
TSUBL = TSUBC * (1 + 0.31 * DF);

% Equation (18)
EXPAP = exp(-0.08 * AP);
DTG = AP + 100 * (1 - EXPAP);

% Compute correction to dTc for local solar time and lat correction
F10 = GEO(1);
DTCLST = DTSUB (F10,GLSTHR,GLAT,ZHT);

% Compute the local exospheric temperature.
TINF = TSUBL + DTG + DTCLST;
TEMP(1) = TINF;

% Equation (9)
TSUBX = 444.3807 + 0.02385 * TINF - 392.8292 * exp(-0.0021357 * TINF);

% Equation (11)
GSUBX = 0.054285714 * (TSUBX - 183.);

% The TC array will be an argument in the call to
% XLOCAL, which evaluates Equation (10) or Equation (13)
TC(1) = TSUBX;
TC(2) = GSUBX;

% A AND GSUBX/A OF Equation (13)
TC(3) = (TINF - TSUBX)/PIOV2;
TC(4) = GSUBX/TC(3);

% Equation (5)
Z1 = 90;
Z2 = min(SAT(3),105);
AL = log(Z2/Z1);
N = floor(AL/R1) + 1;
ZR = exp(AL/N);
AMBAR1 = XAMBAR(Z1);
TLOC1 = XLOCAL(Z1,TC);
ZEND = Z1;
SUM2 = 0;
AIN = AMBAR1 * XGRAV(Z1)/TLOC1;

for I = 1:N
    Z = ZEND;
    ZEND = ZR * Z;
    DZ = 0.25 * (ZEND-Z);
    SUM1 = WT(1)*AIN;
    for J = 2:5
        Z = Z + DZ;
        AMBAR2 = XAMBAR(Z);
        TLOC2 = XLOCAL(Z,TC);
        GRAVL = XGRAV(Z);
        AIN = AMBAR2 * GRAVL/TLOC2;        
        SUM1 = SUM1 + WT(J) * AIN;
    end
    SUM2 = SUM2 + DZ * SUM1;
end

FACT1 = 1000/RSTAR;
RHO = 3.46e-6 * AMBAR2 * TLOC1 * exp(-FACT1*SUM2) /AMBAR1 /TLOC2;

% Equation (2)
ANM = AVOGAD * RHO;
AN  = ANM/AMBAR2;

% Equation (3)
FACT2  = ANM/28.960;
ALN(1) = log(FRAC(1)*FACT2);
ALN(4) = log(FRAC(3)*FACT2);
ALN(5) = log(FRAC(4)*FACT2);

% Equation (4)
ALN(2) = log(FACT2 * (1 + FRAC(2)) - AN);
ALN(3) = log(2 * (AN - FACT2));

if (SAT(3) > 105)
else
    TEMP(2) = TLOC2;
    % Put in negligible hydrogen for use in DO-LOOP 13
    ALN(6) = ALN(5) - 25;    
    % Equation (24)  - J70 Seasonal-Latitudinal Variation
    TRASH = (MJD - 36204) / 365.2422;
    CAPPHI = mod(TRASH,1);
    DLRSL = 0.02 * (SAT(3) - 90) * exp(-0.045 * (SAT(3) - 90))...
          * sign_(1,SAT(2)) * sin(TWOPI * CAPPHI+ 1.72)...
          * sin(SAT(2))^2;
    % Equation (23) - Computes the semiannual variation
    if (Z<2000)
        YRDAY = TMOUTD(MJD);
        % Use new semiannual model
        [FZZ,GTZ,DLRSA] = SEMIAN (YRDAY,ZHT,F10B);
        if (FZZ<0)
            DLRSA = 0;
        end
    end
    
    % Sum the delta-log-rhos and apply to the number densities.
    % In CIRA72 the following equation contains an actual sum,
    % namely DLR = AL10 * (DLRGM + DLRSA + DLRSL)
    % However, for Jacchia 70, there is no DLRGM or DLRSA.
    DLR = AL10 * (DLRSL + DLRSA);
    
    for I = 1:6
        ALN(I) = ALN(I) + DLR;
    end
    
    % Compute mass-density and mean-molecular-weight and
    % convert number density logs from natural to common.
    SUMN = 0;
    SUMNM = 0;
    
    for I = 1:6
        AN = exp(ALN(I));
        SUMN = SUMN + AN;
        SUMNM = SUMNM + AN*AMW(I);
        AL10N(I) = ALN(I)/AL10;
    end
    
    RHO = SUMNM/AVOGAD;
    
    % Compute the high altitude exospheric density correction factor
    FEX = 1;
    
    if ((ZHT>=1000)&&(ZHT<1500))
        ZETA   = (ZHT - 1000) * 0.002;
        ZETA2  =  ZETA * ZETA;
        ZETA3  =  ZETA * ZETA2;
        F15C   = CHT(1) + CHT(2)*F10B + CHT(3)*1500 + CHT(4)*F10B*1500;
        F15C_ZETA = (CHT(3) + CHT(4)*F10B) * 500;
        FEX2   = 3 * F15C - F15C_ZETA - 3;
        FEX3   = F15C_ZETA - 2 * F15C + 2;
        FEX    = 1 + FEX2 * ZETA2 + FEX3 * ZETA3;
    end
    
    if (ZHT >= 1500)
        FEX = CHT(1) + CHT(2)*F10B + CHT(3)*ZHT + CHT(4)*F10B*ZHT;
    end
    
    % Apply the exospheric density correction factor.
    RHO = FEX * RHO;
    return
end

% Equation (6)
Z3 = min(SAT(3),500);
AL = log(Z3/Z);
N = floor(AL/R2) + 1;
ZR = exp(AL/N);
SUM2 = 0;
AIN = GRAVL/TLOC2;

for I = 1:N
    Z = ZEND;
    ZEND = ZR * Z;
    DZ = 0.25 * (ZEND - Z);
    SUM1 = WT(1) * AIN;
    for J = 2:5
        Z = Z + DZ;
        TLOC3 = XLOCAL(Z,TC);
        GRAVL = XGRAV(Z);
        AIN = GRAVL/TLOC3;
        SUM1 = SUM1 + WT(J) * AIN;
    end
    SUM2 = SUM2 + DZ * SUM1;
end

Z4 = max(SAT(3),500);
AL = log(Z4/Z);
R = R2;

if (SAT(3) > 500)
    R = R3;
end

N = floor(AL/R) + 1;
ZR = exp(AL/N);
SUM3 = 0;

for I=1:N
    Z = ZEND;
    ZEND = ZR * Z;
    DZ = 0.25 * (ZEND - Z);
    SUM1 = WT(1) * AIN;
    for J = 2:5
        Z = Z + DZ;
        TLOC4 = XLOCAL(Z,TC);
        GRAVL = XGRAV(Z);
        AIN = GRAVL/TLOC4;
        SUM1 = SUM1 + WT(J) * AIN;
    end
    SUM3 = SUM3 + DZ * SUM1;
end

if (SAT(3) > 500)
    T500 = TLOC3;
    TEMP(2) = TLOC4;
    ALTR = log(TLOC4/TLOC2);
    FACT2 = FACT1 * (SUM2 + SUM3);
    HSIGN = -1;
else
    T500 = TLOC4;
    TEMP(2) = TLOC3;
    ALTR = log(TLOC3/TLOC2);
    FACT2 = FACT1 * SUM2;
    HSIGN = 1;
end

for I = 1:5
    ALN(I) = ALN(I) - (1 + ALPHA(I)) * ALTR - FACT2 * AMW(I);
end

% Equation (7) - Note that in CIRA72, AL10T5 = log10(T500)
AL10T5 = log10(TINF);
ALNH5 = (5.5 * AL10T5 - 39.40) * AL10T5 + 73.13;
ALN(6) = AL10 * (ALNH5 + 6) + HSIGN * (log(TLOC4/TLOC3) + FACT1 * SUM3 * AMW(6));

% Equation (24)  - J70 Seasonal-Latitudinal Variation
TRASH = (MJD - 36204.) / 365.2422;
CAPPHI = mod(TRASH,1);
DLRSL = 0.02 * (SAT(3) - 90) * exp(-0.045 * (SAT(3) - 90))...
      * sign_(1,SAT(2)) * sin(TWOPI * CAPPHI+ 1.72)...
      * sin(SAT(2))^2;

% Equation (23) - Computes the semiannual variation
if (Z<2000)    
    YRDAY = TMOUTD (MJD);
    % Use new semiannual model
    [FZZ,GTZ,DLRSA] = SEMIAN (YRDAY,ZHT,F10B);
    if (FZZ<0)
        DLRSA = 0;
    end
end

% Sum the delta-log-rhos and apply to the number densities.
% In CIRA72 the following equation contains an actual sum,
% namely DLR = AL10 * (DLRGM + DLRSA + DLRSL)
% However, for Jacchia 70, there is no DLRGM or DLRSA.
DLR = AL10 * (DLRSL + DLRSA);

for I = 1:6
    ALN(I) = ALN(I) + DLR;
end

% Compute mass-density and mean-molecular-weight and
% convert number density logs from natural to common.
SUMN = 0;
SUMNM = 0;

for I = 1:6
    AN = exp(ALN(I));
    SUMN = SUMN + AN;
    SUMNM = SUMNM + AN*AMW(I);
    AL10N(I) = ALN(I)/AL10;
end

RHO = SUMNM/AVOGAD;

% Compute the high altitude exospheric density correction factor
FEX = 1;

if ((ZHT>=1000)&&(ZHT<1500))
    ZETA   = (ZHT - 1000) * 0.002;
    ZETA2  =  ZETA * ZETA;
    ZETA3  =  ZETA * ZETA2;
    F15C   = CHT(1) + CHT(2)*F10B + CHT(3)*1500 + CHT(4)*F10B*1500;
    F15C_ZETA = (CHT(3) + CHT(4)*F10B) * 500;
    FEX2   = 3 * F15C - F15C_ZETA - 3;
    FEX3   = F15C_ZETA - 2 * F15C + 2;
    FEX    = 1 + FEX2 * ZETA2 + FEX3 * ZETA3;
end

if (ZHT >= 1500)
    FEX = CHT(1) + CHT(2)*F10B + CHT(3)*ZHT + CHT(4)*F10B*ZHT;
end

% Apply the exospheric density correction factor.
RHO = FEX * RHO;

end

%***********************************************************************
function XAMBAR = XAMBAR(Z)
% Evaluates Equation (1)

C = [28.15204,-8.5586e-2,+1.2840e-4,-1.0056e-5,-1.0210e-5,+1.5044e-6,+9.9826e-8];
DZ = Z - 100;
AMB = C(7);

for I = 1:6
    J = 7-I;
    AMB = DZ * AMB + C(J);
end

XAMBAR = AMB;

end

%***********************************************************************
function XGRAV2 = XGRAV(Z)
% Evaluates Equation (8)
XGRAV2 = 9.80665/(1 + Z/6356.766)^2;

end

%***********************************************************************
function XLOCAL = XLOCAL(Z,TC)
% Evaluates Equation (10) or Equation (13), depending on Z
DZ = Z - 125;
if (DZ > 0)
    XLOCAL = TC(1) + TC(3) * atan(TC(4)*DZ*(1 + 4.5e-6*DZ^2.5));
    return
end
    XLOCAL = ((-9.8204695e-6 * DZ - 7.3039742e-4) * DZ^2 + 1) * DZ * TC(2) + TC(1);
end

%***********************************************************************
function DTC = DTSUB (F10,XLST,XLAT,ZHT)
%
% COMPUTE dTc correction for Jacchia-Bowman model
%
%    Calling Args:
%    ------------
%    F10       = (I)   F10 FLUX
%    XLST      = (I)   LOCAL SOLAR TIME (HOURS 0-23.999)
%    XLAT      = (I)   XLAT = SAT LAT (RAD)
%    ZHT       = (I)   ZHT = HEIGHT (KM)
%    DTC       = (O)   dTc correction
%
B = [-0.457512297e1, -0.512114909e1, -0.693003609e2,...
      0.203716701e3,  0.703316291e3, -0.194349234e4,...
      0.110651308e4, -0.174378996e3,  0.188594601e4,...
     -0.709371517e4,  0.922454523e4, -0.384508073e4,...
     -0.645841789e1,  0.409703319e2, -0.482006560e3,...
      0.181870931e4, -0.237389204e4,  0.996703815e3,...
      0.361416936e2];

C = [-0.155986211e2, -0.512114909e1, -0.693003609e2,...
      0.203716701e3,  0.703316291e3, -0.194349234e4,...
      0.110651308e4, -0.220835117e3,  0.143256989e4,...
     -0.318481844e4,  0.328981513e4, -0.135332119e4,...
      0.199956489e2, -0.127093998e2,  0.212825156e2,...
     -0.275555432e1,  0.110234982e2,  0.148881951e3,...
     -0.751640284e3,  0.637876542e3,  0.127093998e2,...
     -0.212825156e2,  0.275555432e1];

DTC = 0;
tx  = XLST/24;
ycs = cos(XLAT);
F   = (F10 - 100)/100;

% calculates dTc
if (ZHT>=120 && ZHT<=200)
    H = (ZHT - 200)/50;
    DTC200 = C(17)            + C(18)*tx*ycs      + C(19)*tx^2*ycs...
           + C(20)*tx^3*ycs   + C(21)*F*ycs       + C(22)*tx*F*ycs...
           + C(23)*tx^2*F*ycs;
    sum = C(1) + B(2)*F + C(3)*tx*F     + C(4)*tx^2*F...
        + C(5)*tx^3*F    + C(6)*tx^4*F    + C(7)*tx^5*F...
        + C(8)*tx*ycs     + C(9)*tx^2*ycs  + C(10)*tx^3*ycs...
        + C(11)*tx^4*ycs + C(12)*tx^5*ycs + C(13)*ycs...
        + C(14)*F*ycs     + C(15)*tx*F*ycs  + C(16)*tx^2*F*ycs;
    DTC200DZ = sum;
    CC  = 3*DTC200 - DTC200DZ;
    DD  = DTC200 - CC;
    ZP  = (ZHT-120)/80;
    DTC = CC*ZP*ZP + DD*ZP*ZP*ZP;
end

if (ZHT>200 && ZHT<=240)
    H = (ZHT - 200)/50;
    sum = C(1)*H + B(2)*F*H + C(3)*tx*F*H     + C(4)*tx^2*F*H...
        + C(5)*tx^3*F*H    + C(6)*tx^4*F*H    + C(7)*tx^5*F*H...
        + C(8)*tx*ycs*H     + C(9)*tx^2*ycs*H  + C(10)*tx^3*ycs*H...
        + C(11)*tx^4*ycs*H + C(12)*tx^5*ycs*H + C(13)*ycs*H...
        + C(14)*F*ycs*H     + C(15)*tx*F*ycs*H  + C(16)*tx^2*F*ycs*H...
        + C(17)             + C(18)*tx*ycs      + C(19)*tx^2*ycs...
        + C(20)*tx^3*ycs   + C(21)*F*ycs       + C(22)*tx*F*ycs...
        + C(23)*tx^2*F*ycs;
    DTC = sum;
end

if (ZHT>240 && ZHT<=300)
    H = (40)/50;
    sum = C(1)*H + B(2)*F*H + C(3)*tx*F*H     + C(4)*tx^2*F*H...
        + C(5)*tx^3*F*H    + C(6)*tx^4*F*H    + C(7)*tx^5*F*H...
        + C(8)*tx*ycs*H     + C(9)*tx^2*ycs*H  + C(10)*tx^3*ycs*H...
        + C(11)*tx^4*ycs*H + C(12)*tx^5*ycs*H + C(13)*ycs*H...
        + C(14)*F*ycs*H     + C(15)*tx*F*ycs*H  + C(16)*tx^2*F*ycs*H...
        + C(17)             + C(18)*tx*ycs      + C(19)*tx^2*ycs...
        + C(20)*tx^3*ycs   + C(21)*F*ycs       + C(22)*tx*F*ycs...
        + C(23)*tx^2*F*ycs;
    AA = sum;
    BB = C(1) + B(2)*F  + C(3)*tx*F       + C(4)*tx^2*F...
       + C(5)*tx^3*F    + C(6)*tx^4*F    + C(7)*tx^5*F...
       + C(8)*tx*ycs     + C(9)*tx^2*ycs  + C(10)*tx^3*ycs...
       + C(11)*tx^4*ycs + C(12)*tx^5*ycs + C(13)*ycs...
       + C(14)*F*ycs     + C(15)*tx*F*ycs  + C(16)*tx^2*F*ycs;
    H   = 300/100;
    sum = B(1)    + B(2)*F  + B(3)*tx*F         + B(4)*tx^2*F...
        + B(5)*tx^3*F      + B(6)*tx^4*F      + B(7)*tx^5*F...
        + B(8)*tx*ycs       + B(9)*tx^2*ycs    + B(10)*tx^3*ycs...
        + B(11)*tx^4*ycs   + B(12)*tx^5*ycs   + B(13)*H*ycs...
        + B(14)*tx*H*ycs    + B(15)*tx^2*H*ycs + B(16)*tx^3*H*ycs...
        + B(17)*tx^4*H*ycs + B(18)*tx^5*H*ycs + B(19)*ycs;
    DTC300 = sum;
    sum = B(13)*ycs...
         + B(14)*tx*ycs    + B(15)*tx^2*ycs + B(16)*tx^3*ycs...
         + B(17)*tx^4*ycs + B(18)*tx^5*ycs;
    DTC300DZ = sum;
    CC = 3*DTC300 - DTC300DZ - 3*AA - 2*BB;
    DD = DTC300 - AA - BB - CC;
    ZP  = (ZHT-240)/60;
    DTC = AA + BB*ZP + CC*ZP*ZP + DD*ZP*ZP*ZP;
end

if (ZHT>300 && ZHT<=600)
    H   = ZHT/100;
    sum = B(1)    + B(2)*F  + B(3)*tx*F         + B(4)*tx^2*F...
        + B(5)*tx^3*F      + B(6)*tx^4*F      + B(7)*tx^5*F...
        + B(8)*tx*ycs       + B(9)*tx^2*ycs    + B(10)*tx^3*ycs...
        + B(11)*tx^4*ycs   + B(12)*tx^5*ycs   + B(13)*H*ycs...
        + B(14)*tx*H*ycs    + B(15)*tx^2*H*ycs + B(16)*tx^3*H*ycs...
        + B(17)*tx^4*H*ycs + B(18)*tx^5*H*ycs + B(19)*ycs;
        DTC = sum;
end

if (ZHT>600 && ZHT<=800)
    ZP = (ZHT - 600)/100;
    HP = 600/100;
    AA  = B(1)    + B(2)*F  + B(3)*tx*F         + B(4)*tx^2*F...
        + B(5)*tx^3*F      + B(6)*tx^4*F      + B(7)*tx^5*F...
        + B(8)*tx*ycs       + B(9)*tx^2*ycs    + B(10)*tx^3*ycs...
        + B(11)*tx^4*ycs   + B(12)*tx^5*ycs   + B(13)*HP*ycs...
        + B(14)*tx*HP*ycs   + B(15)*tx^2*HP*ycs+ B(16)*tx^3*HP*ycs...
        + B(17)*tx^4*HP*ycs + B(18)*tx^5*HP*ycs + B(19)*ycs;
    BB  = B(13)*ycs...
        + B(14)*tx*ycs    + B(15)*tx^2*ycs + B(16)*tx^3*ycs...
        + B(17)*tx^4*ycs + B(18)*tx^5*ycs;
    CC  = -(3*AA+4*BB)/4.;
    DD  = (AA+BB)/4;
    DTC = AA + BB*ZP + CC*ZP*ZP + DD*ZP*ZP*ZP;
end

end

%******************************************************************
function [FZZ,GTZ,DRLOG] = SEMIAN(DAY,HT,F10BAR)
%
% COMPUTE SEMIANNUAL VARIATION (DELTA LOG RHO)
% INPUT DAY, HEIGHT, F10BAR
%       025.  650.   150.
% OUTPUT FUNCTIONS FZ, GT, AND DEL LOG RHO VALUE
%
% DAY     (I)   DAY OF YEAR
% HT      (I)   HEIGHT (KM)
% F10BAR  (I)   AVE 81-DAY CENTERED F10
% FZZ     (O)   SEMIANNUAL AMPLITUDE
% GTZ     (O)   SEMIANNUAL PHASE FUNCTION
% DRLOG   (O)   DELTA LOG RHO

TWOPI = 2*pi;

% FZ GLOBAL MODEL VALUES
% 1978-2004 FIT:
FZM = [0.111613,-0.159000e-2, 0.126190e-1,...
      -0.100064e-1,-0.237509e-4, 0.260759e-4];

% GT GLOBAL MODEL VALUES
% 1978-2004 FIT:
GTM = [-0.833646,-0.265450, 0.467603,-0.299906,...
       -0.105451,-0.165537e-1,-0.380037e-1,-0.150991e-1,...
       -0.541280e-1, 0.119554e-1, 0.437544e-2,-0.369016e-2,...
        0.206763e-2,-0.142888e-2,-0.867124e-5, 0.189032e-4,...
        0.156988e-3, 0.491286e-3,-0.391484e-4,-0.126854e-4,...
        0.134078e-4,-0.614176e-5, 0.343423e-5];

F10B = F10BAR;
FB2  = F10BAR^2;
HTZ = HT/1000;

FZZ = FZM(1) + FZM(2)*F10B  + FZM(3)*F10B*HTZ + FZM(4)*F10B*HTZ^2 + ...
      FZM(5)*F10B^2*HTZ + FZM(6)*F10B^2*HTZ^2;

TAU   = (DAY-1)/365;
SIN1P = sin(TWOPI*TAU);
COS1P = cos(TWOPI*TAU);
SIN2P = sin(2*TWOPI*TAU);
COS2P = cos(2*TWOPI*TAU);
SIN3P = sin(3*TWOPI*TAU);
COS3P = cos(3*TWOPI*TAU);
SIN4P = sin(4*TWOPI*TAU);
COS4P = cos(4*TWOPI*TAU);

GTZ = GTM(1) + GTM(2)*SIN1P + GTM(3)*COS1P + GTM(4)*SIN2P + GTM(5)*COS2P...
    + GTM(6)*SIN3P + GTM(7)*COS3P + GTM(8)*SIN4P + GTM(9)*COS4P...
    + GTM(10)*F10B + GTM(11)*F10B*SIN1P + GTM(12)*F10B*COS1P...
    + GTM(13)*F10B*SIN2P + GTM(14)*F10B*COS2P + GTM(15)*F10B*SIN3P...
    + GTM(16)*F10B*COS3P + GTM(17)*F10B*SIN4P + GTM(18)*F10B*COS4P...
    + GTM(19)*FB2  + GTM(20)*FB2 *SIN1P + GTM(21)*FB2 *COS1P...
    + GTM(22)*FB2 *SIN2P + GTM(23)*FB2 *COS2P;

if (FZZ<1e-6)
    FZZ = 1e-6;
end

DRLOG = FZZ*GTZ;

end

%******************************************************************
function doy = TMOUTD(MJD)

[year,month,day,hr,min,sec] = invjday(MJD+2400000.5);
doy = finddays(year,month,day,hr,min,sec);

end

