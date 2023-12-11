function [IPop, GMM, IPop_interp] = MASTERIPop(VAR)

% Read MASTER tables
T_LMRO = table2cell(readtable('IPOP_200km_2000km.xlsx','Sheet','LMRO','Range','A2:F36180'));
T_N = table2cell(readtable('IPOP_200km_2000km.xlsx','Sheet','N','Range','A2:E181'));
T_NU = table2cell(readtable('IPOP_200km_2000km.xlsx','Sheet','NU','Range','A2:F36180'));

% Altitudes according to the data download from MASTER
h_min = VAR.h_min;         % [km]
h_max = VAR.h_max;         % [km]
N_shells_master = 181;     % for 10km bin size
R_master = linspace(h_min,h_max,N_shells_master);
alt_meas = 205:10:1995;
re = VAR.re;
R0 = (re+R_master)*1000;
V = 4/3*pi*(diff(R0.^3));   % volume of the shells [m^3]
Vkm = V*1e-9;               % volume of the shells [km^3]

% Number of diameters from MASTER discretisation
N_diam = 200;

% Vector of diameter
d_vect = transpose([T_LMRO{1:N_diam,2}]);
dn_vect = transpose([T_NU{1:N_diam,2}]);

% Create matrices
pop_lmro = zeros(N_shells_master-1,N_diam+1);   % add one column for the altitude
pop_n = zeros(N_shells_master-1,2);             % add one column for the altitude
pop_nu = zeros(N_shells_master-1,N_diam+1);     % add one column for the altitude
k = 1;

for i = 1:N_shells_master-1
    for j = 2:N_diam+1
        pop_lmro(i,j) = [T_LMRO{k,end}];
        pop_nu(i,j) = [T_NU{k,end}];
        k = k+1;
    end
    k = k+1;
end
pop_n(:,2) = [T_N{:,end}];

% From spatial density to number of objects
pop_lmro = pop_lmro .* Vkm';
pop_n = pop_n .* Vkm';
pop_nu = pop_nu .* Vkm';

% Add the altitude column corresponding to MASTER discretisation
pop_lmro(:,1) = alt_meas;
pop_n(:,1) = alt_meas;
pop_nu(:,1) = alt_meas;

% Construct the initial population

% Parameters
N_species = 5; 
ind_S = 3;  ind_D = 4;  ind_N = 5;  ind_B = 6;   ind_NU = 7;
bin_pop = zeros(N_shells_master-1, N_species+2);
bin_pop(:,1) = transpose(alt_meas);

% Subdivide the lmro and untracked debris
for i = 1:N_shells_master-1

    % Bin altitudes
    bin_pop(i,1) = R_master(i);
    bin_pop(i,2) = R_master(i+1);

    for j = 1:N_diam
        if d_vect(j) <= dB_cutoff
            bin_pop(i,ind_S) = bin_pop(i,ind_S) + pop_lmro(i,j+1);
            % bin_pop(i,ind_N) = bin_pop(i,ind_N) + pop_n(i,j+1);                               % to consider rocket bodies also the big fragments
        else
            bin_pop(i,ind_B) = bin_pop(i,ind_B) + pop_lmro(i,j+1);
            % bin_pop(i,ind_S) = bin_pop(i,ind_S) + pop_lmro(i,j+1) + pop_n(i,j+1);             % to consider rocket bodies also the big fragments

        end
        % bin_pop(i,ind_N) = bin_pop(i,ind_N) + pop_n(i,j+1);                                   % to eliminate if we consider as B the big fragments

        if dn_vect(j) >= dN_cutoff(1) && dn_vect(j) <= dN_cutoff(2)
            bin_pop(i,ind_NU) = bin_pop(i,ind_NU) + pop_nu(i,j+1);
        end
    end
end

% Debris
bin_pop(:,ind_N) = pop_n(:,end);

% Distinguish between S and D
D_perc = 0.8;
bin_pop(:,ind_D) = bin_pop(:,ind_S) * D_perc;
bin_pop(:,ind_S) = bin_pop(:,ind_S) - bin_pop(:,ind_D);
bin_pop = round(bin_pop);

% Discretisation
N_shell = (h_max-h_min)/h_bin;
R = linspace(h_min, h_max, N_shell + 1);
IPop = zeros(N_shell, N_species+2);
IPop(:,1) = R(1:end-1);
IPop(:,2) = R(2:end);

for k = 3:N_species+2
    if h_bin == 10
        IPop(:,k) = bin_pop(:,k);
    else
        j = 1;
        for i = 1:size(bin_pop,1)
            if alt_meas(i) <= R(j+1) && alt_meas(i) < h_max
                IPop(j,k) = IPop(j,k) + bin_pop(i,k);
            else
                j = j + 1;
            end
        end
    end
end

% Interpolation
h_interp = 1;
N_interp = (h_max-h_min) / h_interp;
R_interp = linspace(h_min,h_max,N_interp + 1);
IPop_interp = zeros(N_interp, N_species+2);
IPop_interp(:,1) = R_interp(1:end-1);
IPop_interp(:,2) = R_interp(2:end);
IPop_interp(:,3:end) = transpose(spline(IPop(:,2),IPop(:,3:end)',IPop_interp(:,2)));

% GMM
options = statset('MaxIter',500);
for j = 3:N_species+2
    for i = 1:size(IPop,1)
        if i == 1
            GMM.gm{j-2} = IPop(i,2)*ones(IPop(i,j),1);
        else
            GMM.gm{j-2} = [GMM.gm{j-2}; IPop(i,2)*ones(IPop(i,j),1)];
        end
    end
    GMM.gmdist{j-2} = fitgmdist(GMM.gm{j-2},2,'Options',options,'RegularizationValue',0.01);
    GMM.gmsigma{j-2} = GMM.gmdist{j-2}.Sigma;
    GMM.gmmu{j-2} = GMM.gmdist{j-2}.mu;
    GMM.gmwt{j-2} = GMM.gmdist{j-2}.ComponentProportion;
end

% Determine the output
switch spec
    case '3'
        IPop = IPop(:,1:ind_N);
    case '4B'
        IPop = IPop(:,1:ind_B);
    case '5N'
        IPop = IPop(:,1:ind_NU);
    case '5S'
        IPop = 0;
        fprintf('Not implemented\n')
end

end


