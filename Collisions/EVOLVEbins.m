function [nums, isCatastrophic, binOut] = EVOLVEbins(m1,m2,r1,r2,dv,binC,binE,binW,LBdiam,RBflag,sto)
% INPUT
%   m1,m2: mass of objects 1,2 [kg]
%   r1,r2: radii of objects 1,2 [m]
%   dv:    collision vel [km/s]
%   binC:  bin center for mass binning (can be empty: see Bin Definitions below)
%   binE:  bin edges for mass binning (can be empty: see Bin Definitions)
%   binW:  bin widths for mass binning (can be empty: see Bin Definitions)
%   LBdiam: Lower bound of Characteristic Length
%   RBflag: for area to mass ratio (func_Am.m), 1: RB; 0: not RB (default) (optional)
%   sto:   stochastic flag (default: 1) (optional) (0 for deterministic, not implemented yet) 

% OUTPUT
%   nums:  number of objects per bin [N x 1]
%   isCatastrophic:  the collision was catastrophic [1] or not [0]
%   binOut:   bin edges or bin centers depending on input (see below)

% BIN DEFINITIONS
%   1) binC defined  (binE and binW empty: [])
%   2) binC and binW defined  (binE empty)
%   3) binE defined  (binC and binW empty)
%   4) binC and LB defined

% EXAMPLES
%   m1 = 250; m2 = 100; r1 = 0.2616; r2 = 0.1928; dv = 10; binC = [0.6, 2, 10, 50, 100, 200, 300]
%   [nums, isCatastrophic, binOut] = EVOLVEbins(m1,m2,r1,r2,dv,binC,[],[],0.01,RBflag,sto) % Center given
%   [nums, isCatastrophic, binOut] = EVOLVEbins(m1,m2,r1,r2,dv,[],binC,[],0.01,RBflag,sto); % Edges given
%   [nums, isCatastrophic, binOut] = EVOLVEbins(m1,m2,r1,r2,dv,binC,[],0.1,0.1,RBflag,sto); % Center, width given (gap exists)

% Derived from frag_col_SBM_vec.m 
% uses func_Am
% More info on above at the end of this function

SS = 10;  % super sampling ratio (allows for fractional numbers per bin if > 1)

if nargin == 10
    RBflag = 0;
    sto = 1;
elseif nargin == 11
    sto = 1;
else
    error('Wrong number of arguments provided (10 or 11 needed but %i provided)',nargin)
end

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


catastrophRatio = (m2*(dv*1000)^2)/(2*m1*1000);  % J/g = kg*(km/s)^2 / g
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
