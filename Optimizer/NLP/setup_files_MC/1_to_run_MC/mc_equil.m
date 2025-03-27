function [] = mc_equil(seed)

addpath(genpath('../../'));  % add all subfolders of parent folder to path
rng(seed);

if seed <20
    fname = 'mocat3_equilDV_high_cap_1.mat';
    nyears = 100;      
end

fprintf('Seed %i \t %s \t %i yrs\n', seed, fname, nyears);

curfolder = pwd;

fprintf('%s is being loaded in %s \n' , fname, curfolder);
load(fname,'R02','N_shell', 'diameter', 'mass', 'xopt', 'alpha', 'alpha_active', 'P', 'Dt', 'Area');
  
cfgMC = setupMCeqDynamic(seed, nyears, R02, N_shell, diameter, mass, xopt, alpha, alpha_active, P, Dt, Area);

disp('Equilibrium finding...  starting...');
fprintf('Initial Population:  %i sats\n', size(cfgMC.mat_sats,1));
fprintf('Launches per year: %i\n', size(cfgMC.repeatLaunches,1));
disp('Starting main_mc_LNTiadc ...');
main_mc_LNTiadc2(cfgMC);
