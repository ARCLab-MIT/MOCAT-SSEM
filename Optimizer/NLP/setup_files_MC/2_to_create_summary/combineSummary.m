function combineSummary(varargin)

clear
result_dir = pwd; % Store the Results directory path

addpath(genpath('../../'));
getidx;

which Fast_MC2SSEM_population_binned2
which MC2SSEM_population_dist_binned2

% fnames = dir('*.mat');
fnames = dir(fullfile(result_dir, '**', '*.mat')); % '**' for recursive, fullfile for path
fnames = fnames(~[fnames.isdir]);

seeds = cellfun(@(x) round(str2double( x{end-2} )), cellfun(@(x) strsplit(x,'_'),{fnames.name},'UniformOutput',false));

parts = cell2mat(cellfun(@(x) round(str2double(x{end-1})), cellfun(@(x) strsplit(x,{'_','.'}),{fnames.name},'UniformOutput',false),'UniformOutput',false));

% maxparts = splitapply(@max,parts,seeds+1);
maxparts = splitapply(@max,parts,findgroups(seeds));


seededges = [0:20:20]; 
strDescription = {'mocat3_equilDV_high_cap_1'};

% RADIUS BIN VERSION
NmassEdges = []; 
NradiusEdges = [0.05, 0.06,0.07,0.08, 0.10, 0.15, 0.3, 0.45, 0.6, 0.8, inf]; % <<<<<<<<<<<<<<<<<<<<< % 10 species

clear N_MC S_MC D_MC B_MC S N D B pmean pvar pmedian pmeans pvars pmedians pprc  param_mean param_median param_var param_prc
tic
for ind = 2:numel(seededges)
    % curfiles = fnames(seeds >= seededges(ind-1) & seeds < seededges(ind));
    curfiles = {};
    for iseed = seededges(ind-1):seededges(ind)-1
        curseeds = seeds == iseed;
        curparts = parts(curseeds);
        % curmaxpart = maxparts(iseed+1);
        curmaxpart = max(curparts);
	% curmaxpart = 67; % TEMPORARY
        curfiles = [curfiles, fnames(iseed == seeds & curmaxpart == parts)];
    end
    
    foldername  = ['radius' sprintf('_%0.1f',NradiusEdges(2:end-1))];
    mkdir(foldername);
    disp(foldername);
    fi5s = [];
    
    for ii = 1:numel(curfiles)
        clear paramSSEM tmax tmaxx S D N pmeans pmedians pvars pprcs ms sats_info
        load(curfiles{ii}.name)
        fprintf('%s \t x%s | %i/%i seededge | %i/%i files \n',curfiles{ii}.name, strDescription{ind-1}, ind, numel(seededges), ii, numel(curfiles));

        % mass bin edges for Debris [kg]
        paramSSEM.NmassEdges = NmassEdges;
        paramSSEM.NradiusEdges = NradiusEdges;
        % data is from saving every X time steps (parameter used for MC)
        paramSSEM.saveInterval = MCconfig.saveMSnTimesteps;  % 60; should also be in MCconfig.saveMSnTimesteps

        tmax = size(matsatsperN,2);
        tmaxx = tmax;
        tmaxx = 21; % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        fprintf('tmax: %i; tmaxx: %i \n',tmax, tmaxx);
        
        sats_info = cell(1,4);         % contains info for SSEM binning
        S = nan(tmaxx,numel(paramSSEM.R02)-1);
        D = S; B = S;
        % N = nan(tmax,numel(paramSSEM.R02)-1, numel(paramSSEM.NmassEdges)-1);
        N = nan(tmaxx,numel(paramSSEM.R02)-1, numel(paramSSEM.NradiusEdges)-1);
        fprintf('S: %i x %i, D: %i x %i, B: %i x %i, N: %i x %i x %i \n', size(S),size(D),size(B),size(N));
        
        pmeans = nan(paramSSEM.N_shell * (numel(paramSSEM.NradiusEdges)+1), 3 , tmaxx);
        pmedians = pmeans;
        pvars = pmeans;
        pprcs = {};
        
        
        % FOR mc_main_LNTiadc version after 12/6/2023 (separated matsatsperN)
         % merge first using all parts:
            clear MSN
            MSN = {};
            ss = strsplit(curfiles{ii}.name,{'_','.'});
            maxpart = str2num(ss{end-1});
            for pi = 1:maxpart
                loadfilepart = join(ss(1:end-2),'_');
                loadfilepart = [loadfilepart{1}, '_', num2str(pi), '.mat'];
                % disp(loadfilepart)
                load(loadfilepart,'matsatsperN');
                fprintf('successfully loaded matsatsperN in: %s \n', loadfilepart);
                validMSinds = find(~cellfun(@isempty, matsatsperN));
                MSN(validMSinds) = matsatsperN(validMSinds);
            end
            % disp(MSN)
            matsatsperN = MSN;  clear MSN;
        
        disp(size(matsatsperN))
        
             tmax = size(matsatsperN,2);
             fprintf('\t tmax %i\n',tmax)
        
        for tind = 1:tmax
            ms = matsatsperN{tind};
            if isempty(ms)
                S(tind,:) = 0;
                D(tind,:) = 0;
                N(tind,:,:) = 0;
                B(tind,:) = 0;
                % pmeans(:,:,tind) = 0;
                % pmedians(:,:,tind) = 0;
                % pvars(:,:,tind) = 0;
                % % pprcs(1:size(cur_param_prc,1),1:size(cur_param_prc,2),tind) = 0;
                % fprintf('\t tind %i, tmax %i, empty\n',tind,tmax)
                continue;
            else
                % fprintf('\t tind %i, tmax %i, FILLED \n',tind,tmax)
                % CREATE SAT_INFO

                % XXYYZ = p1_objectclass * 1000 + p1_constel * 10 + p1_cont;
                % XXYYZ = mat_sats(:,idx_objectclass) * 1000 + mat_sats(:,idx_constel) * 10 + mat_sats(:,idx_controlled);
                % matsatsperN{1} = single([mat_sats(:,[idx_a, idx_ecco, idx_mass, idx_radius, idx_bstar, idx_ID]),XXYYZ]);
                % OLD: matsatsperN{1} = mat_sats(:,[idx_a, idx_ecco, idx_mass, idx_radius, idx_bstar,  idx_objectclass]);

                sats_info{1} = floor(ms(:,7)./1000); % int8(objclassint_store);
                sats_info{2} = ms(:,1);             % SMA
                sats_info{3} = mod(ms(:,7),10);     % controlled_store
                sats_info{4} = ms(:,3);             % mass
                sats_info{5} = ms(:,4);             % radius

                % CONVERT TO POPSSEM (S_MC, ETC) via Fast_MC2SSEM_population
                %     paramSSEM.NmassEdges = [0,inf];  % default
                % disp('0');
                [S(tind,:), D(tind,:), N(tind,:,:), B(tind,:)]  = Fast_MC2SSEM_population_binned2(sats_info,paramSSEM);

                % CALCULATE PARAM_MEAN/MEDIAN/VAR here (MC2SSEM_population_dist_binned) - EXCLUDES RBs
                [~, cur_param_mean, cur_param_var, cur_param_median, cur_param_prc] = MC2SSEM_population_dist_binned2(ms,paramSSEM);            
                pmeans(:,:,tind) = cur_param_mean;
                pmedians(:,:,tind) = cur_param_median;
                pvars(:,:,tind) = cur_param_var;
                % disp('00');
                % disp(size(pprcs));
                % disp(size(cur_param_prc));
                pprcs(1:size(cur_param_prc,1),1:size(cur_param_prc,2),tind) = cur_param_prc;  % cell
                % disp('1');
            end
        end
        % disp('2');
        S_MC(1:size(S,1),1:size(S,2),ii) = S;
        D_MC(1:size(D,1),1:size(D,2),ii) = D;
        B_MC(1:size(B,1),1:size(B,2),ii) = B;
        N_MC(1:size(N,1),1:size(N,2),1:size(N,3),ii) = N;
        % disp('3');
        param_mean(1:size(pmeans,1),1:size(pmeans,2),1:size(pmeans,3),ii) = pmeans;
        param_median(1:size(pmedians,1),1:size(pmedians,2),1:size(pmedians,3),ii) = pmedians;
        param_var(1:size(pvars,1),1:size(pvars,2),1:size(pvars,3),ii) = pvars;
        param_prc(1:size(pprcs,1),1:size(pprcs,2),1:size(pprcs,3),ii) = pprcs;
        % disp('4');
        % Optional: visualize
        % vizSDNB(S_MC(:,:,1), D_MC(:,:,1), squeeze(N_MC(:,:,1,1))); % for first entry
        if exist('frag_info5') == 1
            fi5s(:,ii) = cellfun(@(x) size(x,1) ,frag_info5);
        end
        toc
    end
        savefilename = sprintf('summary_%s.mat',strDescription{ind-1});
        % savefilename = sprintf('summary_uniform_binned_lamX%0.1f.mat',multipliers(ind-1));

        save([foldername '/' savefilename],'S_MC','D_MC','N_MC','B_MC','param','paramSSEM','MCconfig','param_mean','param_median','param_var','param_prc','fi5s');
        disp([foldername '/' savefilename]);
end
toc

end
