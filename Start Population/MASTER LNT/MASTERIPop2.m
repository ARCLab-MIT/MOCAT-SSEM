function [species_master_quantities, species_master_quantities2] = MASTERIPop2(VAR, deb_species_name, MAX_DIAMETER_CUTOFF, DENSITY)

% Takes scen_properties class. Uses the information associated with the
% shells for that class to apportion master data to the right shells. Also
% accesses the assembled list of debris species to apportion master objects
% to the right species.
% inputs:
%   VAR = a fully constructed scen_properties object with species
%   information.
%   deb_species_name = string or string array corresponding the name of the species in
%   VAR.species_cell that contains the debris population to which MASTER
%   object fluxes will be assigned.
% outputs:
%   species_master_quantities = N_shell x N_species in
%   VAR.species_cell.(deb_species_name) that can be used to augment
%   starting population based on MASTER Launch and Mission Related Objects.
%   species_master_quantities2 = N_shell x N_species in
%   VAR.species_cell.(deb_species_name) that can be used to augment
%   starting population based on MASTER Nu (natural?) objects.

if (~exist('MAX_DIAMETER_CUTOFF', 'var'))
    MAX_DIAMETER_CUTOFF = 0.10; %m
end
if (~exist('DENSITY', 'var'))
    DENSITY = 2710; %km/m^3 Used to estimate mass from diameter assuming Al.
end

% Read MASTER tables
T_LMRO = readtable('IPOP_200km_2000km.xlsx','Sheet','LMRO','ReadVariableNames',true, 'VariableNamingRule', 'preserve' );
T_LMRO = T_LMRO(~any(ismissing(T_LMRO),2),:); %Remove summary stats and blank rows

T_N = readtable('IPOP_200km_2000km.xlsx','Sheet','N','ReadVariableNames',true, 'VariableNamingRule', 'preserve');
T_N = T_N(~any(ismissing(T_N),2),:); %Remove summary stats and blank rows

T_NU = readtable('IPOP_200km_2000km.xlsx','Sheet','NU','ReadVariableNames',true, 'VariableNamingRule', 'preserve');
T_NU = T_NU(~any(ismissing(T_NU),2),:); %Remove summary stats and blank rows

% Count number of shells in file
ind = ~isnan(T_LMRO.Altitude);
N_shells_master = length(unique(T_LMRO.Altitude(ind)));

% Number of diameters from MASTER discretisation
N_diam = sum(T_LMRO.Altitude == T_LMRO.Altitude(1));

% Altitudes according to the data download from MASTER
h_min = min(T_LMRO.Altitude);         % [km]
h_max = max(T_LMRO.Altitude);         % [km]
R_master = linspace(h_min,h_max,N_shells_master+1);
alt_meas = unique(T_LMRO.Altitude(ind));
re = VAR.re;
R0 = (re+R_master)*1000;
V = 4/3*pi*(diff(R0.^3));   % volume of the shells [m^3]
Vkm = (V*1e-9);               % volume of the shells [km^3]

% Vector of diameter
d_vect = transpose([T_LMRO{1:N_diam,2}]);
dn_vect = transpose([T_NU{1:N_diam,2}]);

%%
% Create matrices
pop_lmro = zeros(N_shells_master,N_diam+1);   % add one column for the altitude
pop_n = zeros(N_shells_master,2);             % add one column for the altitude
pop_nu = zeros(N_shells_master,N_diam+1);     % add one column for the altitude

for i = 1:length(alt_meas)
    alt = alt_meas(i);
    T_LMRO_total_col = T_LMRO(T_LMRO.Altitude == alt,["Diameter", "Total"]);
    T_LMRO_total_row = rows2vars(T_LMRO_total_col, 'VariableNamesSource', ...
                                 'Diameter', 'VariableNamingRule', 'preserve');
    T_LMRO_total_row = renamevars(T_LMRO_total_row,'OriginalVariableNames','Altitude');
    T_LMRO_total_row.Altitude = alt;
    pop_lmro(i,:) = table2array(T_LMRO_total_row);

    T_NU_total_col = T_NU(T_NU.Altitude == alt,["Diameter", "Total"]);
    T_NU_total_row = rows2vars(T_NU_total_col, 'VariableNamesSource', ...
                                 'Diameter', 'VariableNamingRule', 'preserve');
    T_NU_total_row = renamevars(T_NU_total_row,'OriginalVariableNames','Altitude');
    T_NU_total_row.Altitude = alt;
    pop_nu(i,:) = table2array(T_NU_total_row);
end

pop_n(:,2) = [T_N{:,"Total"}];

% From spatial density to number of objects
pop_lmro = pop_lmro .* Vkm';
pop_n = pop_n .* Vkm';
pop_nu = pop_nu .* Vkm';

% Add the altitude column corresponding to MASTER discretisation
pop_lmro(:,1) = alt_meas;
pop_n(:,1) = alt_meas;
pop_nu(:,1) = alt_meas;

pop_lmro = array2table(pop_lmro, 'VariableNames', ["Altitude", num2cell(d_vect)] );
pop_nu = array2table(pop_nu, 'VariableNames', ["Altitude", num2cell(dn_vect)] );

% Construct the initial population
pop_lmro_model_shells = zeros(VAR.N_shell, length(d_vect));

% Find objects in model shells based on VAR.
for i = 1:length(VAR.R02)-1
    shell_bot = VAR.R02(i);
    shell_top = VAR.R02(i+1);
    %disp(["Bottom: ", shell_bot, "Top: ", shell_top])
    shell_rows = pop_lmro((pop_lmro.Altitude >= shell_bot) & (pop_lmro.Altitude < shell_top),:);
    shell_total = sum(shell_rows{:,2:end},1);
    shell_total = array2table(shell_total, 'VariableNames', shell_rows.Properties.VariableNames(2:end));
    pop_lmro_model_shells(i,:) = table2array(shell_total);

    shell_rows2 = pop_nu((pop_nu.Altitude >= shell_bot) & (pop_nu.Altitude < shell_top),:);
    shell_total2 = sum(shell_rows2{:,2:end},1);
    shell_total2 = array2table(shell_total2, 'VariableNames', shell_rows2.Properties.VariableNames(2:end));
    pop_nu_model_shells(i,:) = table2array(shell_total2);
end

% Turn diameter into mass values based on spherical aluminum assumption
v_d_vect = (d_vect./2).^3 * (4/3) * pi;
m_d_vect = v_d_vect .* DENSITY;

v_dn_vect = (dn_vect./2).^3 * (4/3) * pi;
m_dn_vect = v_dn_vect .* DENSITY;

MAX_MASS_CUTOFF = (MAX_DIAMETER_CUTOFF/2)^3 * (4/3) * pi * DENSITY; %m

% Get the relevant species from VAR
species_cell = species.empty;
for name_i = 1:length(deb_species_name)
    getfield(VAR.species_cell, deb_species_name(name_i));
    species_cell = [species_cell getfield(VAR.species_cell, deb_species_name(name_i))];
end

species_master_quantities = zeros(VAR.N_shell, length(species_cell));
species_master_quantities2 = zeros(VAR.N_shell, length(species_cell));
species_names = strings(3,1);

for i = 1:length(species_cell)
    species_names(i,1) = species_cell(i).species_properties.sym_name;
    mass_lb = species_cell(i).species_properties.mass_lb;

    % Only want to get objects from master up to size chosen by user
    mass_ub = min(species_cell(i).species_properties.mass_ub, MAX_MASS_CUTOFF);
    
    in_mass_band = pop_lmro_model_shells(:, (mass_lb<=m_d_vect) & (m_d_vect < mass_ub));
    in_mass_band_total = sum(in_mass_band, 2);
    species_master_quantities(:,i) = in_mass_band_total;

    in_mass_band2 = pop_nu_model_shells(:, (mass_lb<=m_dn_vect) & (m_dn_vect < mass_ub));
    in_mass_band_total2 = sum(in_mass_band2, 2);
    species_master_quantities2(:,i) = in_mass_band_total2;
end

species_master_quantities = array2table(species_master_quantities, 'VariableNames', species_names);
species_master_quantities2 = array2table(species_master_quantities2, 'VariableNames', species_names);


