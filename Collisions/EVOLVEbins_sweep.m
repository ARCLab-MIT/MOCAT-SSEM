% sweep thru EVOLVEbins

% Example (from EVOLVEbins.m): 
%   m1 = 250; m2 = 100; r1 = 0.2616; r2 = 0.1928; dv = 10; binC = [0.6, 2, 10, 50, 100, 200, 300]
%   [nums, isCatastrophic, binOut] = EVOLVEbins(m1,m2,r1,r2,dv,binC,[],[],0.01,RBflag,sto) % Center given


% sweep masses and DVs
m1s = 10:10:50;  % kg
dvs = 1:12;  % km/s 

m2 = 1; % kg
r1 = 3; % m
r2 = 2; % m - may not matter tho

binc = [0.5, 1, 5 : 5 : 50];  % mass bin centers for output (kg)
RBflag = 0;   % assume not RB
sto = 0;

sweepN = nan(numel(m1s),numel(dvs));
for mind = 1:numel(m1s)
    for dvind = 1:numel(dvs)
        m1 = m1s(mind);
        dv = dvs(dvind);
        [nums, edges, iscat] = EVOLVEbins(m1,m2,r1,r2,dv,binc,[],[],0.01,RBflag,sto);
%                              EVOLVEbins(m1,m2,r1,r2,dv,binC,[],[],0.01,RBflag,sto) % Center given
        sweepN(mind,dvind) = sum(nums);
%         fprintf('mind: %i, dvind: %i, cat: %i \n',mind,dvind,iscat);
    end
end

% Sweep result (total number of frags)
figure(2); clf;
surf(dvs,m1s,sweepN)
xlabel('DV [km/s]')
ylabel('Object 1 mass [kg]')
zlabel('Total number of fragments');
title('Total number of fragments')

%% One example of binned frags
figure(3); clf; 
X = categorical(cellstr(num2str(binc')));
% X = categorical({'a','b','c','d'});
X = reordercats(X,cellstr(num2str(binc')));
% X = reordercats(X, {'a','b','c','d'});
b = bar(X, nums);
xlabel('Bin centers (kg)');
ylabel('Debris count in bin');
title('EVOLVE Bin distribution for', ...
    sprintf('m1: %0.1f; m2: %0.2f kg; r1: %0.2f; r2; %0.2f m; dv: %0.1f km/s',...
    m1,m2,r1,r2,dv));

%% SURF plot: num of objects for [m1 mass] x [kg bins]
% sweep masses
m1s = 10:50:500;  % kg
dv = 10;  % km/s 

m2 = 1; % kg
r1 = 3; % m
r2 = 1; % m - may not matter tho

% binc = [0.5, 1, 2 : 2 : 30];  % mass bin centers for output (kg)
binc = logspace(-1,1.3,10)/2;

RBflag = 0;   % assume not RB
sto = 0;

sweepN = nan(numel(m1s),numel(binc));
for mind = 1:numel(m1s)
    m1 = m1s(mind);
    [nums, iscat, binout] = EVOLVEbins(m1,m2,r1,r2,dv,binc,[],[],0.01,RBflag,sto);  % Center given
    sweepN(mind,:) = nums;
    fprintf('mind: %0.1f, iscat: %i \n',mind, iscat);
end

% Sweep result (total number of frags)
figure(2); clf;
% surf(binc,m1s,sweepN)
surf(sweepN);
a = gca;
a.XTickLabel = num2str(binc','%0.2f');
a.YTickLabel = num2str(m1s');
xlabel('Debris mass bins [kg]','interpreter','tex','fontsize',15);
ylabel('Parent 1 mass [kg]','interpreter','tex','fontsize',15);
zlabel('Number of fragments per bin w_k','interpreter','tex','fontsize',15);
% title('Total number of fragments')

% print("binningSweep",'-depsc2','-vector');

%% SURF plot #2: For a given P1 vs P2 distribution, how does it look different for various binC's?
% aka, go DEEP into EVOLVEbins.m  to deconstruct within (center version only):

m1 = 200;  % kg
m2 = 1; % kg
r1 = 3; % m
r2 = 2; % m - may not matter tho
RBflag = 0;   % assume not RB
sto = 0;

dv = 10;  % km/s 

LB = 0.01;

% function [nums, isCatastrophic, binOut] = EVOLVEbins(m1,m2,r1,r2,dv,binC,binE,binW,LBdiam,RBflag,sto)

SS = 10;  % super sampling ratio (allows for fractional numbers per bin if > 1)

% DO COLLISIONS HERE?
objclass = 5;  % RB for func_Am.m
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
    error('0 debris generated');
%     nums = zeros(numel(binEd)-1,1);
%     binOut = [];
%     return
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
dss = d_pdf(randi(numel(d_pdf),ceil(numSS),1)); % Limit number of fragments to be equal to 'num'


% calculate mass of objects [LB, 1 m] by d > A > Am > m
A = 0.556945*dss.^(2.0047077);      % calculate area; Eq 2.72
Am = func_Am(dss, objclass);     % use Am conversion of larger object <<<<
m = A./Am;

%% histogram of raw m's  (alt 1 and 2 in mass_binning.pptx)
figure(10); clf;
subplot(311) % length
% h = histogram(dss,20);  a = gca; a.YScale = 'log'; x = xlim; xlim([0,x(2)])
% h.BinCounts = h.BinCounts ./ SS;
% xlabel('Length of debris (m)'); ylabel('Count')
% subplot(312) % mass

h = histogram(m,20); a = gca; a.YScale = 'log'; x = xlim; xlim([0,x(2)])
h.BinCounts = h.BinCounts ./ SS;
xlabel('Mass of debris (m)'); ylabel('Count')

% subplot(313); % various bins
subplot(312); % various bins
binC = [0.5, 1, 5 : 5 : 20];  % mass bin centers for output (kg)
% LBm = binC(1) - (binC(2) - binC(1))/2;  % extrapolate for lowest edge
LBm = 0; % lowest edge is always 0 kg
UBm = binC(end) + (binC(end) - binC(end-1))/2; % for highest edge
binEd = [LBm, mean([ binC(1:end-1); binC(2:end)]), UBm];
nums = histcounts(m,binEd) / SS;  % desample here; bin via histcounts
histogram('BinCounts', nums, 'BinEdges', binEd,'edgealpha',0.5,'FaceAlpha',0.15); a = gca; a.YScale = 'log';
xlabel('Mass of debris (m)'); ylabel('Count')
hold on;
xlim([0,x(2)])


subplot(313); % various bins
binC = [0.5, 2, 7, 20];  % mass bin centers for output (kg)
% LBm = binC(1) - (binC(2) - binC(1))/2;  % extrapolate for lowest edge
UBm = binC(end) + (binC(end) - binC(end-1))/2; % for highest edge
binEd = [LBm, mean([ binC(1:end-1); binC(2:end)]), UBm];
nums = histcounts(m,binEd) / SS;  % desample here; bin via histcounts
histogram('BinCounts', nums, 'BinEdges', binEd,'edgealpha',0.5,'FaceAlpha',0.1, ...
    'facecolor','r'); a = gca; a.YScale = 'log';
xlabel('Mass of debris (m)'); ylabel('Count')
xlim([0,x(2)])


% histogram2('BinCounts', nums, 'XBinEdges', binEd, 'YBinEdges', [1,2]); a = gca; a.YScale = 'log';
% xlabel('Mass of debris (m)'); ylabel('Count')

% histogram2('XBinEdges',-1:1,'YBinEdges',-2:2,'BinCounts',[1 2 3 4; 5 6 7 8])


% figure;
% E = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.1, 2.7, 3.6, 5.1];
% C = [10, 80, 132, 142, 138, 112, 75, 65, 52, 65, 62, 32, 35];
% histogram('BinCounts', C, 'BinEdges', E)

% histogram2('XBinEdges',Xedges,'YBinEdges',Yedges,'BinCounts',counts) manually specifies the bin counts. histogram2 plots the specified bin counts and does not do any data binning.

%% histogram of raw m's  (alt 3)

figure(10); clf;
subplot(311) % length
h = histogram(m,100); a = gca; a.YScale = 'log'; 
%x = xlim; xlim([0,x(2)])
xlim([0,18]); x = xlim;
h.BinCounts = h.BinCounts ./ SS;
xlabel('Mass of debris (kg)'); ylabel('Count'); ylim([0.05,1e5])

subplot(312); % various bins
binC = [0.5, 1, 5 : 5 : 15];  % mass bin centers for output (kg)
% LBm = binC(1) - (binC(2) - binC(1))/2;  % extrapolate for lowest edge
LBm = 0; % lowest edge is always 0 kg
UBm = binC(end) + (binC(end) - binC(end-1))/2; % for highest edge
binEd = [LBm, mean([ binC(1:end-1); binC(2:end)]), UBm];
nums = histcounts(m,binEd) / SS;  % desample here; bin via histcounts
n = [zeros(1,numel(nums)); nums];
n = n(:);
b = [binC - 0.1; binC + 0.1]; bb = [0; b(:)];

histogram('BinCounts', n', 'BinEdges', bb','edgealpha',0.5,'FaceAlpha',0.5, ...
    'facecolor','g'); a = gca; a.YScale = 'log'; xlabel('Mass bins (kg)'); ylabel('Count')
hold on;
xlim([0,x(2)]); ylim([0.05,1e5])
a = gca; a.XTick = binC;

subplot(313); % various bins
binC = [0.5, 3, 7, 15];  % mass bin centers for output (kg)
% LBm = binC(1) - (binC(2) - binC(1))/2;  % extrapolate for lowest edge
UBm = binC(end) + (binC(end) - binC(end-1))/2; % for highest edge
binEd = [LBm, mean([ binC(1:end-1); binC(2:end)]), UBm];
nums = histcounts(m,binEd) / SS;  % desample here; bin via histcounts
n = [zeros(1,numel(nums)); nums];
n = n(:);
b = [binC - 0.1; binC + 0.1]; bb = [0; b(:)];
histogram('BinCounts', n', 'BinEdges', bb','edgealpha',0.5,'FaceAlpha',0.5, ...
    'facecolor','r'); a = gca; a.YScale = 'log';
xlabel('Mass bins (kg)'); ylabel('Count')
xlim([0,x(2)]); ylim([0.05,1e5])
a = gca; a.XTick = binC;
