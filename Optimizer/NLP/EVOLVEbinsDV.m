function [nums, isCatastrophic, binOut, altNums] = EVOLVEbinsDV(m1,m2,r1,r2,dv,binC,binE,binW,LBdiam,RBflag,R02)

% INPUT
%   m1,m2: mass of objects 1,2 [kg]
%   r1,r2: radii of objects 1,2 [m]
%   dv:    collision vel [km/s]
%   binC:  bin center for mass binning (can be empty: see Bin Definitions below)
%   binE:  bin edges for mass binning (can be empty: see Bin Definitions)
%   binW:  bin widths for mass binning (can be empty: see Bin Definitions)
%   LBdiam: Lower bound of Characteristic Length
%   RBflag: for area to mass ratio (func_Am.m), 1: RB; 0: not RB (default) (optional)
%   R02:   altitude bin edges (km)

% OUTPUT
%   nums:  number of objects per bin [N x 1]
%   isCatastrophic:  the collision was catastrophic [1] or not [0]
%   binOut:   bin edges or bin centers depending on input (see below)
%   altNums: 2*nShells-1 x 1   of distribution of number of debris into each altitude shell
%                           [-(nShells-1):(nShells-1)]

% BIN DEFINITIONS
%   1) binC defined  (binE and binW empty: [])
%   2) binC and binW defined  (binE empty)
%   3) binE defined  (binC and binW empty)
%   4) binC and LB defined

% EXAMPLES
%   m1 = 250; m2 = 100; r1 = 0.2616; r2 = 0.1928; dv = 10; 
%   binC = [0.6, 2, 10, 50, 100, 200, 300]; R02 = 200:50:2000;
%   [nums, isCatastrophic, binOut] = EVOLVEbins(m1,m2,r1,r2,dv,binC,[],[],0.01,RBflag,R02) % Center given
%   [nums, isCatastrophic, binOut] = EVOLVEbins(m1,m2,r1,r2,dv,[],binC,[],0.01,RBflag,R02); % Edges given
%   [nums, isCatastrophic, binOut] = EVOLVEbins(m1,m2,r1,r2,dv,binC,[],0.1,0.1,RBflag,R02); % Center, width given (gap exists)

% Derived from frag_col_SBM_vec.m 
% uses func_Am
% More info on above at the end of this function

SS = 20;  % super sampling ratio (allows for fractional numbers per bin if > 1)

% if nargin == 10
%     RBflag = 0;
%     sto = 1;
% elseif nargin == 11
%     sto = 1;
% else
%     error('Wrong number of arguments provided (10 or 11 needed but %i provided)',nargin)
% end

% Define and validate bin edges
if ~isempty(binC) && isempty(binE) && isempty(binW)         % option 1: bin center given
    % LB = 0.1;       %10 cm lower bound; L_c [m]
    LBm = binC(1) - (binC(2) - binC(1))/2;  % extrapolate for lowest edge
    UBm = binC(end) + (binC(end) - binC(end-1))/2; % for highest edge
    binEd = [LBm, mean([ binC(1:end-1); binC(2:end)]), UBm];
elseif ~isempty(binC) && ~isempty(binW) && isempty(binE)    % option 2: bin width given; bin gaps exist
    if any(diff(binC) < binW)
        errinds = find(diff(binC) < binW);
        error('Overlapping bin edges between bin centered at %0.1f and %0.1f',...
            binC(errinds(1)), binC(errinds(1)+1));
    end
    binEd1 = binC - binW/2;
    binEd2 = binC + binW/2;
    binEd = sort([binEd1, binEd2]);
elseif ~isempty(binE) && isempty(binC) && isempty(binW)     % option 3: bin edges given
    binEd = sort(binE);
else
    error('Wrong setup for bins given (binC empty: %i; binE empty: %i; binW empty: %i)',...
        isempty(binC),isempty(binE),isempty(binW));
end

LB = LBdiam;

if RBflag == 0
    objclass = 5;  % RB for func_Am.m
else
    objclass = 0;  % not RB
end

% Ensure p1_mass > p2_mass
if m2 > m1
    tm = m2;
    m2 = m1;
    m1 = tm;
    tr = r2;
    r2 = r1;
    r1 = tr;
end


catastrophRatio = (m2*(dv*1000)^2)/(2*m1*1000);  % J/g = kg*(m/s)^2 / g
% If the specific energy is < 40 J/g: non-catastrophic collision
if catastrophRatio<40
    M = m2*dv^2;           % correction from ODQN [kg*km^2/s^2]
    isCatastrophic = 0;
else    % catastrophic collision
    M = m1 + m2;
    isCatastrophic = 1;
end

% num = floor(0.1 * M ^ (0.75) * LB ^(-1.71) - 0.1 * M ^ (0.75) * min([1,2*r1]) ^(-1.71));
num = (0.1 * M ^ (0.75) * LB ^(-1.71)) - (0.1 * M ^ (0.75) * min([1,2*r1]) ^(-1.71));
numSS = SS * num;
%fprintf('Collision created %i objects\n', num);

if numSS == 0  % check if 0 (e.g. if LB = r1)
    nums = zeros(numel(binEd)-1,1);
    binOut = [];
    return
end

% Create PDF of power law dist, then sample 'num' selections
% only up to 1m, then randomly sample larger objects as quoted above

dd_edges = logspace(log10(LB),log10(min([1,2*r1])),500); % log space, up to either 1 m or diameter of larger satellite

log10_dd = log10(dd_edges); %log10 of diameter edge bins
dd_means = 10.^(log10_dd(1:end-1)+diff(log10_dd)/2); %mean value of each diameter edge bin (in log scale, not linear scale, since bins are generated in a logarithmic scale)

nddcdf = 0.1 * M^(0.75) * dd_edges.^(-1.71);  % CUMULATIVE distribution for collision  (eq 2.68)
ndd = max(0,-diff(nddcdf)); % diff to get PDF count for the bins (dd_edges); if negative, set to zero
floor_ndd = floor(ndd); %floor of PDF count for each bin
rand_sampling = rand(size(ndd)); %0 to 1, random number for stochastic sampling of fragment diameters
add_sampling = rand_sampling>(1-(ndd-floor_ndd)); %0 if random number is lower than 1-(decimal part of ndd), 1 if random number is above 1-(decimal part)
d_pdf = repelem(dd_means,floor_ndd+add_sampling)';   % PDF of debris objects between LB and 1 m

% d = d_pdf(randperm(numel(d_pdf))); % Do not limit number of fragments to be equal to 'num'
try
    dss = d_pdf(randi(numel(d_pdf),ceil(numSS),1)); % Limit number of fragments to be equal to 'num'
catch % When the probability distribution breaks due to small objects colliding.
    % disp("DSS errored for m1: " + string(m1) + " m2: " + string(m2) + ...
    %      " r1: " + r1 + " r2: " + string(r2) + " dv: " + string(dv))
    % disp("binc: " + string(binC))
    % disp("bine: ")
    % disp(string(binE)) 
    % disp("binw: " + string(binW))
    % disp("LBdiam: " + string(LBdiam) + " RBflag: " + string(RBflag) ...
    %     + " sto: " + string(sto))
    dss = 0;
end

% calculate mass of objects [LB, 1 m] by d > A > Am > m
A = 0.556945*dss.^(2.0047077);      % calculate area; Eq 2.72
Am = func_Am(dss, objclass);     % use Am conversion of larger object <<<<
m = A./Am;

% do the binning via histcounts
nums = histcounts(m,binEd) / SS;  % desample here

% Option 2: remove numbers corresponding to bin gaps
if ~isempty(binW) 
    fprintf('Removing %i objects from bin gaps (out of %i total)\n',...
        sum(nums(2:2:end)), sum(nums));
    nums(2:2:end) = [];
end

% define binOut -- changes depending on option
binOut = [];  % option 2
if ~isempty(binC) && isempty(binE) && isempty(binW)         % option 1: bin center given; output = edges
    binOut = binEd;
elseif ~isempty(binE) && isempty(binC) && isempty(binW)     % option 3: bin edges given; output = centers
    binOut = binE(1:end-1) + diff(binE)/2;
end


%
%
% Assign delta-v to spherically random directions
    % test: DV difference between shells
    % MU = 398600.4418;  RE = 6378.1;  % km
    % figure;plot(200:50:2000 ,sqrt( MU ./ (RE+ (200:50:2000))),'-x'); grid on;
    % figure;plot(25 + [200:50:1950], diff(sqrt( MU ./ (RE+(200:50:2000)))*1000),'-x'); ylabel('m/s'); grid on;
    
    dAlt = median(diff(R02));  % diff between shells (km)
    nShell = numel(diff(R02));
    
    % find difference in orbital velocity for shells
    MU = 398600.4418;  RE = 6378.1;  % km
    % use equal spacing in DV space for binning to altitude bins
    dDV = abs(median(diff(sqrt( MU ./ (RE+(200:50:2000)))*1000))); % m/s original
    % dDV = abs(median(diff(sqrt( MU ./ (RE+(200:50:900)))*1000)));   % m/s

    dv = func_dv(Am,'col')/1000;       % km/s
    u=rand(length(dv),1)*2-1;
    theta=rand(length(dv),1)*2*pi;
    
    v = sqrt(1 - u.^2);
    p = [v.*cos(theta) v.*sin(theta) u];
    dv_vec = p .* dv; 
    hc = histcounts(dv_vec(:),[-(nShell-0.5):nShell-0.5]*dDV/1000); % merge directions to oversample by SS * 3
    altNums = hc ./ SS ./ 3;  % may not sum to nums bc of <R02(1) or >R02(end)
    

% references
%   Johnson, N. L., et al. "NASA's New Breakup Model of EVOLVE 4.0" 2001
%   Klinkrad, "Space Debris: Models and Risk Analysis" 2006
%   MASTER-8-Final-Report
%   https://fragmentation.esoc.esa.int/home/modelling


% Key idea
% 1) determine if catastrophic collision (> 40 J/g)
% 2) N(d) = 0.1 * m_e ^ 0.75 * d^-1.71
%       if catastrophic: m_e = total mass
%       if non-catastrophic: m_e = m_p*v_i^2
% 3) Sample bottom-up until "total mass" as defined in the quote below is
%    reached. Create one remnant object (non-catastrophic), or several
%    (catastrophic) large fragements with remnant mass

end

function z=func_dv(Am,s)

% reference: SPACE DEBRIS, Model and Risk Analysis
% eq(3.42)

switch s
    case 'exp'
        mu = 0.2*log10(Am) + 1.85;  % explosion
    case 'col'
        mu = 0.9*log10(Am) + 2.9;   % collision
    otherwise
        warning('exp/col not specified; using explosion parameter')
end


sigma = 0.4;

N = mu+sigma*randn(size(mu));
% N = normrnd(mu,sigma,1,1);

z = 10.^N;        % m/s
end

function out = func_Am(d,ObjClass)
% NASA's new breakup model of evolve 4.0
% Author: N.L.Johnson, P.H.Krisko, J.-C. Liou, P.D.Anz-Meador
% Eq(6) pg 1381
% Distribution function for spacecraft fragments with Lc(=d)>11cm 

numObj = numel(d);      
logds = log10(d);        % d in meters
amsms = nan(numObj,5);  % store alpha,mu1,sig1,mu2,sig2

% if strcmpi(ObjClass,'Rocket Body') ||...
%    strcmpi(ObjClass,'Rocket Mission Related Object') ||...
%    strcmpi(ObjClass,'Rocket Fragmentation Debris') ||...
%    strcmpi(ObjClass,'Rocket Debris')
if ObjClass>4.5 && ObjClass<8.5 %Rocket-body related
    for ind = 1:numObj
        logd = logds(ind);
        % alpha
        if logd <= -1.4
            alpha = 1;
        elseif (-1.4 < logd) && (logd < 0)
            alpha = 1 - 0.3571 * (logd + 1.4);
        else % >= 0
            alpha = 0.5;
        end
    
        % mu1
        if logd <= -0.5
            mu1 = -0.45;
        elseif (-0.5 < logd) && (logd < 0)
            mu1 = -0.45 - 0.9 * (logd + 0.5);
        else % >= 0
            mu1 = -0.9;
        end
    
        % sigma1
        sigma1 = 0.55;
    
        % mu2
        mu2 = -0.9;
    
        % sigma2
        if logd <= -1.0
            sigma2 = 0.28;
        elseif (-1 < logd) && (logd < 0.1)
            sigma2 = 0.28 - 0.1636 * (logd + 1);
        else % >= 0.1
            sigma2 = 0.1;
        end
        amsms(ind,:) = [alpha, mu1, sigma1, mu2, sigma2];
    end
else        % not rocket body
    for ind = 1:numObj
        logd = logds(ind);
        % alpha
        if logd <= -1.95
            alpha = 0;
        elseif (-1.95 < logd) && (logd < 0.55)
            alpha = 0.3 + 0.4 * (logd + 1.2);
        else % >= 0.55
            alpha = 1;
        end
    
        % mu1
        if logd <= -1.1
            mu1 = -0.6;
        elseif (-1.1 < logd) && (logd < 0)
            mu1 = -0.6 - 0.318 * (logd + 1.1);
        else % >= 0
            mu1 = -0.95;
        end
    
        % sigma1
        if logd <= -1.3
            sigma1 = 0.1;
        elseif (-1.3 < logd) && (logd < -0.3)
            sigma1 = 0.1 + 0.2 * (logd + 1.3);
        else % >= -0.3
            sigma1 = 0.3;
        end
    
        % mu2
        if logd <= -0.7
            mu2 = -1.2;
        elseif (-0.7 < logd) && (logd < -0.1)
            mu2 = -1.2 - 1.333 * (logd + 0.7);
        else % >= -0.1
            mu2 = -2.0;
        end
    
        % sigma2
        if logd <= -0.5
            sigma2 = 0.5;
        elseif (-0.5 < logd) && (logd < -0.3)
            sigma2 = 0.5 - (logd + 0.5);
        else % >= -0.3
            sigma2 = 0.3;
        end
        amsms(ind,:) = [alpha, mu1, sigma1, mu2, sigma2];
    end
end

% N1 = mu1+sigma1*randn;
N1 = amsms(:,2) + amsms(:,3) .* randn(numObj,1);
% N2 = mu2+sigma2*randn;
N2 = amsms(:,4) + amsms(:,5) .* randn(numObj,1);

% out = (10.^(alpha * N1 + (1-alpha) * N2)); % eq(3.40)
out = (10.^(amsms(:,1) .* N1 + (1-amsms(:,1)) .* N2));

% Unit, seems to be m^2/kg per Fig 6.

end
