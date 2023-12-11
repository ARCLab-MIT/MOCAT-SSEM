function dens_lowvar = dens_prediction(n_years, start_year, start_month, start_day)

    % dens_prediction - Generate jb2008 density field for n_years from
    % chosen start period
    % 2020
    %
    % This code is licensed under the GNU General Public License version 3.
    %
    % Author: Peng Mun Siew
    % Massachusetts Institute of Technology, Dept. of Aeronautics and
    % Astronautics
    % email: siewpengmun@yahoo.com
    % Jun 2022; Last revision: 22-Nov-2022
    
    
    n_years = 3;
    start_year = 2020;
    start_month = 3;
    start_day = 15;
    
    % clearvars
    % clc
    %n_years = 200;
    %TODO: Bounds checking from files
    
    %% LOAD KERNELS, GRAVITY MODEL, EARTH ORIENTATION PARAMETERS AND SGP4
    % Load SPICE
    fileName = matlab.desktop.editor.getActiveFilename;
    [filepath,name,ext] = fileparts(fileName);
    
    spicePath = fullfile(filepath, '../','mice'); 
    addpath(fullfile(spicePath,'src','mice'));
    addpath(fullfile(spicePath,'lib'));
    
    % progress bar
    cd ../../../
    progressPath = fullfile(pwd, 'Dependencies', 'MatlabProgressBar');
    addpath(progressPath)
    cd(filepath)

    kernelpath  = fullfile( '.\Data','kernel.txt');
    loadSPICE(kernelpath);
    addpath('JB2008');
    
    % Load solar and geomagnetic indices based on CIRA-08 (Moderately active
    % solar cycle)
    load('yearssolar.mat')
    
    %% F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B
    F10mean = yearssolar.F10mean;
    F10std = min([yearssolar.F10max-yearssolar.F10mean,yearssolar.F10mean-yearssolar.F10min],[],2);
    
    F81mean = yearssolar.F81mean;
    F81std = min([yearssolar.F81max-yearssolar.F81mean,yearssolar.F81mean-yearssolar.F81min],[],2);
    
    S10mean = yearssolar.S10mean;
    S10std = min([yearssolar.S10max-yearssolar.S10mean,yearssolar.S10mean-yearssolar.S10min],[],2);
    
    S81mean = yearssolar.S81mean;
    S81std = min([yearssolar.S81max-yearssolar.S81mean,yearssolar.S81mean-yearssolar.S81min],[],2);
    
    M10mean = yearssolar.M10mean;
    M10std = min([yearssolar.M10max-yearssolar.M10mean,yearssolar.M10mean-yearssolar.M10min],[],2);
    
    M81mean = yearssolar.M81mean;
    M81std = min([yearssolar.M81max-yearssolar.M81mean,yearssolar.M81mean-yearssolar.M81min],[],2);
    
    Y10mean = yearssolar.Y10mean;
    Y10std = min([yearssolar.Y10max-yearssolar.Y10mean,yearssolar.Y10mean-yearssolar.Y10min],[],2);
    
    Y81mean = yearssolar.Y81mean;
    Y81std = min([yearssolar.Y81max-yearssolar.Y81mean,yearssolar.Y81mean-yearssolar.Y81min],[],2);
    
    n_rep = ceil(n_years/12);
    F10_monthly = (repmat(F10std,[n_rep,1])./3).*randn(n_rep*12*12,1) + repmat(F10mean,[n_rep,1]);
    F10B_monthly = (repmat(F81std,[n_rep,1])./3).*randn(n_rep*12*12,1) + repmat(F81mean,[n_rep,1]);
    S10_monthly = (repmat(S10std,[n_rep,1])./3).*randn(n_rep*12*12,1) + repmat(S10mean,[n_rep,1]);
    S10B_monthly = (repmat(S81std,[n_rep,1])./3).*randn(n_rep*12*12,1) + repmat(S81mean,[n_rep,1]);
    M10_monthly = (repmat(M10std,[n_rep,1])./3).*randn(n_rep*12*12,1) + repmat(M10mean,[n_rep,1]);
    M10B_monthly = (repmat(M81std,[n_rep,1])./3).*randn(n_rep*12*12,1) + repmat(M81mean,[n_rep,1]);
    Y10_monthly = (repmat(Y10std,[n_rep,1])./3).*randn(n_rep*12*12,1) + repmat(Y10mean,[n_rep,1]);
    Y10B_monthly = (repmat(Y81std,[n_rep,1])./3).*randn(n_rep*12*12,1) + repmat(Y81mean,[n_rep,1]);
    
    DSTDTC = 58; %TODO: What is this?
    
    %% Month vector
    month_vec = [start_month:12,1:start_month-1]'; %TODO: ASK PENG if this needs to change
    mmUTC_flat = repmat(month_vec,[n_rep*12,1]);
    yyUTC = repmat(start_year:start_year+n_rep*12+1,[12,1]);
    yyUTC_flat = reshape(yyUTC,[],1);
    yyUTC_flat = yyUTC_flat(start_month:n_rep*12*12+2);
    
    yyUTC = yyUTC_flat(1);
    mmUTC = mmUTC_flat(1);
    %% JB2008 density
    
    doyUTC = day(datetime(yyUTC, mmUTC, start_day),'dayofyear');
    
    [monthUTC,dayUTC,~,~,~] = days2mdh(yyUTC,doyUTC);
    MJD = Mjday(yyUTC,monthUTC,dayUTC,0,0,0);
    
    [DJMJD0, DATE] = iauCal2jd(yyUTC, mmUTC, 15);
    % TIME = (60*(60*hour+minute)+sec)/86400;
    % UTC = DATE+TIME;
    UTC = DATE;
    
    et  = cspice_str2et(strcat([num2str(jed2date(UTC+2400000.5),'%d %d %d %d %d %.10f') 'UTC']));
    rr_sun = cspice_spkezr('Sun',et,'J2000','NONE', 'Earth');
    rr_sun = rr_sun(1:3,1);
    ra_Sun  = atan2(rr_sun(2), rr_sun(1));
    dec_Sun = atan2(rr_sun(3), sqrt(rr_sun(1)^2+rr_sun(2)^2));
    
    SUN(1)  = ra_Sun;
    SUN(2)  = dec_Sun;
    
    %%  lat, lon in pi
    
    lon_vec = linspace(0,2*pi,91);
    lat_vec = linspace(-0.98*pi/2,0.98*pi/2,46);
    alt = linspace(200,2000,37); %km
    
    [lon2d,lat2d] = meshgrid(lon_vec,lat_vec);
    
    lon = reshape(lon2d,[],1);
    lat = reshape(lat2d,[],1);
    
    dens_table = zeros(length(alt),length(F10_monthly));
    % month_ik = 1;
    
    %TODO add parfor
    for month_ik = progress(1:length(F10_monthly))
        
        % Update MJD
        MJD = Mjday(yyUTC_flat(month_ik),mmUTC_flat(month_ik),15,0,0,0);
        
        % Update SUN position
        et  = cspice_str2et(strcat([num2str(jed2date(MJD+2400000.5),'%d %d %d %d %d %.10f') 'UTC']));
        rr_sun = cspice_spkezr('Sun',et,'J2000','NONE', 'Earth');
        rr_sun = rr_sun(1:3,1);
        ra_Sun  = atan2(rr_sun(2), rr_sun(1));
        dec_Sun = atan2(rr_sun(3), sqrt(rr_sun(1)^2+rr_sun(2)^2));
        SUN(1)  = ra_Sun;
        SUN(2)  = dec_Sun;
        
        % Update space weather indices
        F10 = F10_monthly(month_ik);
        F10B = F10B_monthly(month_ik);
        S10 = S10_monthly(month_ik);
        S10B = S10B_monthly(month_ik);
        M10 = M10_monthly(month_ik);
        M10B = M10B_monthly(month_ik);
        Y10 = Y10_monthly(month_ik);
        Y10B = Y10B_monthly(month_ik);
    
        for height_ik = 1:length(alt)
            rho_jb2 = zeros(length(lon),1);
    
            for ik = 1:length(lon)
                SAT(1) = lon(ik); % 0: 2*pi - lon_geocentric
                SAT(2) = lat(ik); % -pi : pi - lat_geocentric
                SAT(3) = alt(height_ik); % alt(km)
                [~,rho_jb2(ik)] = JB2008(MJD,SUN,SAT,F10,F10B,S10,S10B,M10,M10B,Y10,Y10B,DSTDTC);
            end
            dens_table(height_ik,month_ik) = mean(rho_jb2);
        end
    end
    
    %% Convert data into a structure
    dens_lowvar.dens = dens_table;
    dens_lowvar.year = yyUTC_flat;
    dens_lowvar.month = mmUTC_flat;
    dens_lowvar.alt = alt;
    
    %% Save data for MC_model
    save('dens_lowvar_2000.mat','dens_lowvar')
    
    %%
    close all
    
    for ik = 1:17
       plot(dens_table(5,1+144*(ik-1):144*ik))
       hold on
    end
    grid on
    xlabel('Month')
    ylabel('Density at 300 km')
    title('Variation of density across different solar cycles')
    
    %%
    figure (2)
    for ik = 1:17
       semilogy(alt,dens_table(:,ik*10))
       hold on
    end
    grid on
    hold off
    xlabel('Altitude')
    ylabel('Density')
    title('Variation of density across altitudes for different months')
