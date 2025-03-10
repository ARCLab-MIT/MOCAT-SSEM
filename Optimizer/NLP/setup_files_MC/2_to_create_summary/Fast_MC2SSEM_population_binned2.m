function [S_MC, D_MC, N_MC, RB_MC] = Fast_MC2SSEM_population_binned2(sats_info,param)

% Divides the output of MC runs (sats_info) by species and altitude bin. 
% The structure param contains the parameters that need to be set. (e.g.  param.NmassEdges or  param.NradiusEdges)
%---
% Authors: Peng Mun Siew, MIT 11/4/2022
%          djang - modified to add binning per param.NmassEdges; to be used with
%                  convert2popSSEM, 7/1/2023
%                  Re-added RBs, 7/31
%---

% species: [S,D,N,RB];

idx_active = sats_info{1}==1 & sats_info{3}==1;
idx_derelict = sats_info{1}==1 & sats_info{3}~=1;
idx_debris = sats_info{1}==3 | sats_info{1}==4 | sats_info{1}>=7;
idx_rbodies = sats_info{1}==5 | sats_info{1}==6;

obj_alt = sats_info{2}* param.re - param.re;

% Separate sat alt based on classes
alt_active = obj_alt(idx_active);
alt_derelict = obj_alt(idx_derelict);
alt_debris = obj_alt(idx_debris);
mass_debris = sats_info{4}(idx_debris);
radius_debris = sats_info{5}(idx_debris);
alt_rbodies = obj_alt(idx_rbodies);

[s_count,~] = histcounts(alt_active,param.R02);
[d_count,~] = histcounts(alt_derelict,param.R02);
% [n_count,~] = histcounts(alt_debris,param.R02);
if ~isempty(param.NmassEdges)
    n_count_bin = histcounts2(mass_debris, alt_debris, param.NmassEdges, param.R02);
    % disp(['Debris class binned by MASS (kg): '  num2str(param.NmassEdges)])
elseif ~isempty(param.NradiusEdges)
    n_count_bin = histcounts2(radius_debris, alt_debris, param.NradiusEdges, param.R02);
    % disp(['Debris class binned by RADIUS (m): '  num2str(param.NradiusEdges)])
else
    error('param.NmassEdges and NradiusEdges are both empty!');
end

[b_count,~] = histcounts(alt_rbodies,param.R02);
% su_count = zeros(1,param.N_shell);
% u_count = zeros(1,param.N_shell);

S_MC = s_count' ;
D_MC = d_count' ;
N_MC = n_count_bin' ;
RB_MC = b_count' ;

% if sum(idx_rbodies) > 0
%     warning('rocket bodies detected in sats_info - these are not counted');
% end
end