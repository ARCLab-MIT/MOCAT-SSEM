function [popSSEM, popSSEM_param_mean, popSSEM_param_var, popSSEM_param_median, popSSEM_param_prc] ...
    = MC2SSEM_population_dist_binned2(ms,param)

% Divides the output of the supercloud runs by species and altitude bin.
% The structure param contains the parameters that need to be set.
%---
% Authors: Peng Mun Siew, MIT 11/4/2022
%          djang 07/01/2023  -- for binned N class objects (see Fast_MC2SSEM_population_binned.m)
%                            -- removed support for RBs
%                07/26       -- add other percentiles (2,10,25,75,90,98th)
%                12/15       -- fixed object definition (XXYYZ) to new matsatsperN definition
%---

% ORIGINAL MATSATS DEFINITION
% idx_a = 1; idx_ecco = 2; idx_inclo = 3; idx_nodeo = 4; idx_argpo = 5; idx_mo = 6; idx_bstar = 7; idx_mass = 8; idx_radius = 9;
% idx_error = 10; idx_controlled = 11; idx_a_desired = 12; idx_missionlife = 13; idx_constel = 14; idx_date_created = 15; idx_launch_date = 16;
% idx_r = [17 18 19]; idx_v = [20 21 22]; idx_objectclass = 23; idx_ID = 24;

% Save flag 10 save method in main_mc2.m: 
% matsatsperN{1} = single(ms(:,[idx_a, idx_ecco, idx_mass, idx_radius, idx_bstar,  idx_objectclass, idx_controlled]));
    idx_a = 1;
    idx_objectclass = 6;
    idx_controlled = 7;
    idx_bstar = 5;
    idx_mass = 3;
    idx_radius = 4;

% species: [S,D,N]


% XXYYZ = p1_objectclass * 1000 + p1_constel * 10 + p1_cont;
% XXYYZ = mat_sats(:,idx_objectclass) * 1000 + mat_sats(:,idx_constel) * 10 + mat_sats(:,idx_controlled);
% matsatsperN{1} = single([mat_sats(:,[idx_a, idx_ecco, idx_mass, idx_radius, idx_bstar, idx_ID]),XXYYZ]);
% OLD: matsatsperN{1} = mat_sats(:,[idx_a, idx_ecco, idx_mass, idx_radius, idx_bstar,  idx_objectclass]);
idx_active = floor(ms(:,7)/1000)==1 & rem(ms(:,7),10)==1;  % payload & controlled
idx_derelict = floor(ms(:,7)/1000)==1 & rem(ms(:,7),10)~=1;
idx_debris = floor(ms(:,7)/1000)==3 | floor(ms(:,7)/1000)==4 | floor(ms(:,7)/1000)>=7;

obj_alt = ms(:,idx_a)* param.re - param.re;

% Separate sat alt based on classes
alt_active = obj_alt(idx_active); %load('summary_uniform_binned_x7.mat')
% vizParamStats(param_mean,param_median,param_var,paramSSEM)

alt_derelict = obj_alt(idx_derelict);
alt_debris = obj_alt(idx_debris);

% Separate sat parameters based on classes
param_active = ms(idx_active,[idx_bstar,idx_mass,idx_radius]);  % nObj x 3
param_derelict = ms(idx_derelict,[idx_bstar,idx_mass,idx_radius]);
param_debris = ms(idx_debris,[idx_bstar,idx_mass,idx_radius]);

% create super structure for debris 
    if ~isempty(param.NmassEdges)             % if binned by MASS
        [~,~,n_bin] = histcounts(param_debris(:,2), param.NmassEdges);
        param_debris_bin = {}; alt_debris_bin = {};
        for ind = 1 : numel(param.NmassEdges)-1
            param_debris_bin{ind} = param_debris(n_bin == ind,:);
            alt_debris_bin{ind} = alt_debris(n_bin == ind);
        end
        % disp(['Debris class binned by MASS (kg): '  num2str(param.NmassEdges)])
        
    elseif ~isempty(param.NradiusEdges)       % if binned by RADIUS
        [~,~,n_bin] = histcounts(param_debris(:,3), param.NradiusEdges);
        param_debris_bin = {}; alt_debris_bin = {};
        for ind = 1 : numel(param.NradiusEdges)-1
            param_debris_bin{ind} = param_debris(n_bin == ind,:);
            alt_debris_bin{ind} = alt_debris(n_bin == ind);
        end
        % disp(['Debris class binned by RADIUS (m): '  num2str(param.NradiusEdges)])
    else
        error('NmassEdges and NradiusEdges are both empty!');
    end
 
[s_count,~,s_bin] = histcounts(alt_active,param.R02);
[d_count,~,d_bin] = histcounts(alt_derelict,param.R02);
% [n_count,~,n_bin] = histcounts(alt_debris,param.R02);
% [n_count_bin,~,~,n_Xbin,n_Ybin] = histcounts2(param_debris(:,2), alt_debris,...
%     param.NmassEdges, param.R02);  % [N,Xedges,Yedges,binX,binY] = histcounts2(___)
n_count_bin = {}; n_bin_bin = {};
for ind = 1:numel(param_debris_bin)
    [nc,~,nb] = histcounts(alt_debris_bin{ind},param.R02);
    n_count_bin{ind} = nc;
    n_bin_bin{ind} = nb;
end

param_active_bin_mean = zeros(param.N_shell,size(param_active,2));
param_active_bin_var =  zeros(param.N_shell,size(param_active,2));
param_active_bin_median = zeros(param.N_shell,size(param_active,2));
param_derelict_bin_mean = zeros(param.N_shell,size(param_active,2));
param_derelict_bin_var =  zeros(param.N_shell,size(param_active,2));
param_derelict_bin_median = zeros(param.N_shell,size(param_active,2));
param_debris_bin_mean = zeros(param.N_shell,size(param_active,2));
param_debris_bin_var =  zeros(param.N_shell,size(param_active,2));
param_debris_bin_median = zeros(param.N_shell,size(param_active,2));
param_active_bin_prc = {};
param_derelict_bin_prc = {};
param_debris_bin_prc = {};

% Remove bin that are out of range
param_active(s_bin==0,:) = [];
param_derelict(d_bin==0,:) = [];
% param_debris(n_bin==0,:) = [];
% param_rbodies(b_bin==0,:) = [];
s_bin(s_bin==0) = [];
d_bin(d_bin==0) = [];
% n_bin(n_bin==0) = [];
for ind = 1:numel(n_bin_bin)
    param_debris_bin{ind}(n_bin_bin{ind} == 0,:) = [];
    n_bin_bin{ind}(n_bin_bin{ind} == 0) = [];    
end

prc = @(x) {prctile(x,[5,25,75,95])}; % percentile function for *_bin_prc

for ik = 1:size(param_active,2)  % num of attributes, e.g. [idx_bstar,idx_mass,idx_radius]
    if ~isempty(s_bin)
        param_active_bin_mean(:,ik) = accumarray(s_bin,param_active(:,ik),[param.N_shell 1],@mean);
        param_active_bin_var(:,ik) = accumarray(s_bin,param_active(:,ik),[param.N_shell 1],@var);
        param_active_bin_median(:,ik) = accumarray(s_bin,param_active(:,ik),[param.N_shell 1],@median);
        param_active_bin_prc(:,ik) = accumarray(s_bin,param_active(:,ik),[param.N_shell 1], prc);
    end
    if ~isempty(d_bin)
        param_derelict_bin_mean(:,ik) = accumarray(d_bin,param_derelict(:,ik),[param.N_shell 1],@mean);
        param_derelict_bin_var(:,ik) = accumarray(d_bin,param_derelict(:,ik),[param.N_shell 1],@var);
        param_derelict_bin_median(:,ik) = accumarray(d_bin,param_derelict(:,ik),[param.N_shell 1],@median);
        param_derelict_bin_prc(:,ik) = accumarray(d_bin,param_derelict(:,ik),[param.N_shell 1], prc);
    end
    for ind = 1:numel(n_bin_bin)
        n_bin = n_bin_bin{ind};
        param_debris = param_debris_bin{ind};
        if ~isempty(n_bin)
            param_debris_bin_mean_bin{ind}(:,ik) = accumarray(n_bin,param_debris(:,ik),[param.N_shell 1],@mean);
            param_debris_bin_var_bin{ind}(:,ik) = accumarray(n_bin,param_debris(:,ik),[param.N_shell 1],@var);
            param_debris_bin_median_bin{ind}(:,ik) = accumarray(n_bin,param_debris(:,ik),[param.N_shell 1],@median);
            param_debris_bin_prc_bin{ind}(:,ik) = accumarray(n_bin,param_debris(:,ik),[param.N_shell 1], prc);
        else
            param_debris_bin_mean_bin{ind}(:,ik) = nan(param.N_shell,1);
            param_debris_bin_var_bin{ind}(:,ik) = nan(param.N_shell,1);
            param_debris_bin_median_bin{ind}(:,ik) = nan(param.N_shell,1);
            param_debris_bin_prc_bin{ind}(:,ik) = cell(param.N_shell,1);
        end
    end
end

% popSSEM = [s_count' d_count' n_count' su_count' b_count' u_count' ];
popSSEM = [s_count' d_count' cell2mat(n_count_bin(:))'];

popSSEM_param_mean = [param_active_bin_mean;param_derelict_bin_mean; vertcat(param_debris_bin_mean_bin{:})]; %;param_rbodies_bin_mean];
popSSEM_param_var = [param_active_bin_var;param_derelict_bin_var; vertcat(param_debris_bin_var_bin{:})]; %;param_rbodies_bin_var];
popSSEM_param_median = [param_active_bin_median;param_derelict_bin_median; vertcat(param_debris_bin_median_bin{:})]; %param_rbodies_bin_median];
popSSEM_param_prc = [param_active_bin_prc;param_derelict_bin_prc; vertcat(param_debris_bin_prc_bin{:})];

% clean up stats' zeros into nans
popSSEM_param_mean(popSSEM_param_mean  == 0) = nan;
popSSEM_param_var(popSSEM_param_var  == 0) = nan;
popSSEM_param_median(popSSEM_param_median  == 0) = nan;
% popSSEM_param_prc(popSSEM_param_prc == 0) = nan;
popSSEM_param_prc(cellfun(@isempty, popSSEM_param_prc)) = {nan(1,4)};

end