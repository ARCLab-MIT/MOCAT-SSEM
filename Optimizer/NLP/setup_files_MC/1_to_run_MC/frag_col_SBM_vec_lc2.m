function [debris1, debris2, isCatastrophic] = frag_col_SBM_vec_lc2(ep, p1_in, p2_in, param, LB)
% Collision model following NASA EVOLVE 4.0 standard breakup model (2001)
% with the revision in ODQN "Proper Implementation of the 1998 NASA Breakup Model" (2011)
%p1_in = [p1.mass,p1.radius,p1.r,p1.v,p1.objectclass]
%p2_in = [p2.mass,p2.radius,p2.r,p2.v,p2.objectclass]

% Mar 2024 - non-catastrophic collision creates parent object with large enough mass as
% "large debris" (missing before)

% inputs
%   ep: epoch
%   p1: object 1 (fields used: mass, radius, v, objectclass)
%   p2: object 2
%   param: holds param.max_frag, param.mu, param.req
%   LB : lower bound (meters); aka characteristic length

% references
%   Johnson, N. L., et al. "NASA's New Breakup Model of EVOLVE 4.0" 2001
%   Klinkrad, "Space Debris: Models and Risk Analysis" 2006
%   MASTER-8-Final-Report

% Key idea
% 1) determine if catastrophic collision (> 40 J/g)
% 2) N(d) = 0.1 * m_e ^ 0.75 * d^-1.71
%       if catastrophic: m_e = total mass
%       if non-catastrophic: m_e = m_p*v_i^2
% 3) Sample bottom-up until "total mass" as defined in the quote below is
%    reached. Create one remnant object (non-catastrophic), or several
%    (catastrophic) large fragements with remnant mass

% Quote (ODQN):
% Catastrophic collision fragments are deposited from 1 mm upward along
% that collisional characteristic length distribution until the bin before 1 m. 
% The total mass, M = Mt + Mp , is achieved through deposit of several large 
% fragments on that last bin. 
% Non-catastrophic collision fragments are deposited from 1 mm upward 
% until the total mass, M = Mp * v_i^2  is achieved.  The final fragment 
% is deposited in a single massive fragment reminiscent of a cratered target mass.

% [ed] Note that the SBM distribution is modified by ESA for MASTER-2009 and MASTER-8
% and also by other authors in other papers e.g. Cimmino (2021), etc.
% These are not implemented in this SBM code.

% "Validation" with the ESA's "Validated" C++ code is on-going. 

% LB = 0.1; %10 cm lower bound; L_c

% Ensure p1_mass > p2_mass, or p1_radius>p2_radius if p1_mass=p2_mass
if p1_in(1)<p2_in(1) || (p1_in(1)==p2_in(1) && p1_in(2)<p2_in(2)) %p1_mass < p2_mass || (p1_mass == p2_mass && p1_radius < p2_radius)
    temp1 = p1_in;   temp2 = p2_in;
    p1_in = temp2;   p2_in = temp1;
end

p1_mass = p1_in(1); p1_radius = p1_in(2);
p1_r = p1_in(3:5); p1_v = p1_in(6:8);
p1_objclass = p1_in(9);
p2_mass = p2_in(1); p2_radius = p2_in(2);
p2_r = p2_in(3:5); p2_v = p2_in(6:8);
p2_objclass = p2_in(9);

dv = norm(p1_v - p2_v);                      % km/s
catastrophRatio = (p2_mass*(dv*1000)^2)/(2*p1_mass*1000);  % J/g = kg*m^2/s^2 / g
% plot catastrophic ratio delineation curve? 

% If the specific energy is < 40 J/g: non-catastrophic collision
if catastrophRatio<40
    M = p2_mass*dv^2;           % correction from ODQN [kg*km^2/s^2]
                                % power-law debris until M, then rest is deposited into
                                % one large fragment from larger parent
    isCatastrophic = 0;
else    % catastrophic collision
    M = p1_mass + p2_mass;
    isCatastrophic = 1;
end

num = floor(0.1 * M ^ (0.75) * LB ^(-1.71) - 0.1 * M ^ (0.75) * min([1,2*p1_radius]) ^(-1.71));
% Create PDF of power law dist, then sample 'num' selections
% only up to 1m, then randomly sample larger objects as quoted above
% NOTE: 'num' ISN'T USED -- the same number is recreated in nddcdf below (binned)

% dd_edges = linspace(LB,10,500);   % sampled diameters ; power law valid from LB to 1 m
% dd_edges = logspace(log10(LB),0,200); % log space instead
dd_edges = logspace(log10(LB),log10(min([1,2*p1_radius])),200); % log space, up to either 1 m or diameter of larger satellite
log10_dd = log10(dd_edges); %log10 of diameter edge bins
dd_means = 10.^(log10_dd(1:end-1)+diff(log10_dd)/2); %mean value of each diameter edge bin (in log scale, not linear scale, since bins are generated in a logarithmic scale)

nddcdf = 0.1 * M^(0.75) * dd_edges.^(-1.71);  % CUMULATIVE distribution for collision  (eq 2.68)
ndd = max(0,-diff(nddcdf)); % diff to get PDF count for the bins (dd_edges); if negative, set to zero
floor_ndd = floor(ndd); %floor of PDF count for each bin
rand_sampling = rand(size(ndd)); %0 to 1, random number for stochastic sampling of fragment diameters
add_sampling = rand_sampling>(1-(ndd-floor_ndd)); %0 if random number is lower than 1-(decimal part of ndd), 1 if random number is above 1-(decimal part)
d_pdf = repelem(dd_means,floor_ndd+add_sampling)';   % PDF of debris objects between LB and 1 m

d = d_pdf(randperm(numel(d_pdf))); % Do not limit number of fragments to be equal to 'num'
% d = d_pdf(randi(numel(d_pdf),ceil(num),1)); % Limit number of fragments to be equal to 'num'
% 
% figure(50); clf;
% histogram(d,dd_edges)
% hold on
% plot(dd_means,ndd,'k')
% set(gca,'YScale','log','XScale','log')
% xlabel('d [m]')
% ylabel('Number of fragments [-]')
% title('Number of fragments vs diameter')
% 
%     figure(22); clf; subplot(321); 
%     Lcs = logspace(-2,1,100); loglog(Lcs, 0.1 * M^(0.75) * Lcs.^(-1.71)); grid on; ylim([1,10000]); 
%     hold on; loglog(dd_edges(1:end),nddcdf); xlabel('diam'); ylabel('CDF'); title('theoretical reverse CDF');
%     legend('NASA: N_{cum} = 0.1 M_p^{0.75} L_c^{-1.71}','nddcdf from code');
%     subplot(323); histogram(d_pdf); xlabel('diam'); ylabel('count'); title('PDF to sample from')
%     subplot(325); loglog(dd_means, flip(cumsum(flip(histcounts(repelem(dd_means,round(ndd))',dd_edges)))),'-x'); xlabel('Diam (m)'); ylabel('cumulative count'); title('CDF of above (for shape)');
%     subplot(324); histogram(d); xlabel('diam'); ylabel('count'); title('PDF of sampled d');
%     subplot(326); loglog(dd_edges, 0.1 * M^(0.75) * dd_edges.^(-1.71)); hold on; loglog(dd_means, flip(cumsum(flip(histcounts(d,dd_edges)))),'-x'); xlabel('Diam (m)'); ylabel('cumulative count'); title('CDF of sampled debris sizes');
    
% calculate mass of objects [LB, 1 m] by d > A > Am > m
A = 0.556945*d.^(2.0047077);            % calculate area; Eq 2.72
Am = func_Am(d, p1_objclass);        % use Am conversion of larger object <<<<
m = A./Am;

% create all fragments first (either from p1 or p2)
if sum(m) < M        % if debris mass is less than total mass,
    if isCatastrophic
        % Catastrophic collision fragments are deposited from 1 mm upward along
        % that collisional characteristic length distribution until the bin before 1 m. 
        % The total mass, M = Mt + Mp , is achieved through deposit of several large 
        % fragments on that last bin.       
        largeidx = (d > p2_radius*2 | m > p2_mass) & d < p1_radius*2; %0s or 1s, 1 if larger than smaller satellite and shorter than larger satellite (except for remnant objects that have already been assigned to one of the satellites)
        m_assigned_large = max([0,sum(m(largeidx))]); %mass of fragments that need to belong to larger satellite
        
        if m_assigned_large>p1_mass %if assigned fragments to larger satellite already amount to more than the mass of the satellite
            % Sort by Lc, keep smallest objects until mass adds up <<<<<<<<<<<<<<
            idx_large = find(largeidx); %find indices of fragments assigned to larger satellite
            [~, dord1] = sort(m(idx_large));
            cumsum_m1 = cumsum(m(idx_large(dord1))); %cumulative sum of fragment masses in ascending order of diameter
            lastidx1 = max([0,find(cumsum_m1 < p1_mass, 1, 'last')]); % last valid index of sorted list
            m(idx_large(dord1(lastidx1+1:end))) = [];   %Remove fragments exceeding mass constraint, that cannot belong to the smaller parent either
            d(idx_large(dord1(lastidx1+1:end))) = [];
            A(idx_large(dord1(lastidx1+1:end))) = [];
            Am(idx_large(dord1(lastidx1+1:end))) = [];
            largeidx(idx_large(dord1(lastidx1+1:end))) = [];
            if lastidx1==0
                m_assigned_large = 0; %redefine mass assigned to larger satellite
            else
                m_assigned_large = cumsum_m1(lastidx1); %redefine mass assigned to larger satellite
            end
        end
        mass_max_small = min([p2_mass,m_assigned_large]); %maximum mass to fill up with fragments in smaller satellite
        
        smallidx_temp = find(~largeidx); %fragments that have not been assigned to larger satellite yet
        [~, dord] = sort(m(smallidx_temp)); %sort fragments in ascending order of diameter, those that do not belong to larger satellite
        cumsum_m = cumsum(m(smallidx_temp(dord))); %cumulative sum of fragment masses in ascending order of diameter, those that do not belong to larger satellite
        lastidx_small = find(cumsum_m <= mass_max_small, 1, 'last');  % last valid index of sorted list for cumulative mass under mass_max_small
%         if ~isempty(lastidx_small) %if at least one fragment should be assigned to smaller satellite        
        smallidx = false(size(d)); %preallocate fragments that should be assigned to smaller satellite
        smallidx(smallidx_temp(dord(1:lastidx_small))) = 1;  %0s or 1s, 1 for fragments that belong to smaller satellite
        m_assigned_small = max([0,sum(m(smallidx))]); %mass of fragments that need to belong to smaller satellite

        m_remaining_large = p1_mass-m_assigned_large; %remaining mass in larger satellite
        m_remaining_small = p2_mass-m_assigned_small; %remaining mass in smaller satellite
        m_remaining = [m_remaining_large,m_remaining_small]; %vector with remaining masses associated to each parent satellite
        
        % Assign remnant mass > d > A > Am > dv
        m_remSum = M - sum(m);                  % remnant mass
        remDist = rand(randi([2,8]),1);         % randomly distribute by random number [2,8] <<<<<<<
        m_rem_temp = m_remSum * remDist / sum(remDist); % randomly distribute the remnant mass
        num_rem = numel(m_rem_temp); %number of remnant objects
        
        [m_rem_sort,idx_m_sort] = sort(m_rem_temp,'descend'); % order remnant objects in descending order of mass
        rem_temp_ordered = 1+round(rand(num_rem,1)); %1s or 2s, randomly assign ordered remnant objects to either of the satellites, 1 for larger, 2 for smaller parent        
        for i_rem=1:num_rem %loop for each remnant object
            if m_rem_sort(i_rem)>m_remaining(rem_temp_ordered(i_rem)) %check if mass of remnant object is larger than remaining mass of corresponding parent
                diff_mass = m_rem_sort(i_rem)-m_remaining(rem_temp_ordered(i_rem)); %mass that is left over from remnant object, that could not be assigned to corresponding parent satellite
                m_rem_sort(i_rem) = m_remaining(rem_temp_ordered(i_rem)); %assign remaining mass of corresponding parent to remnant object
                m_remaining(rem_temp_ordered(i_rem)) = 0; %no remaining mass left from corresponding parent satellite
                if i_rem==num_rem %if this is the last remnant object, create one extra remnant object with remaining mass and assign to the other parent
                    m_rem_sort(end+1) = diff_mass; %assign remaining mass to extra remnant object
                    idx_m_sort(end+1) = num_rem+1; %include one extra index in sorted vector of masses
                    rem_temp_ordered(end+1) = ((rem_temp_ordered(i_rem)-1)==0)+1; %1 or 2, assign extra remnant object to the other parent
                else %distribute remaining mass from remmant object amongst other remnant objects
                    m_rem_sort(i_rem+1:end) = m_rem_sort(i_rem+1:end)+diff_mass/(num_rem-i_rem); %equally distribute remaining mass amongst other remnant objects
                    rem_temp_ordered(i_rem+1:end) = ((rem_temp_ordered(i_rem)-1)==0)+1; %1 or 2, assign remaining remnant objects to the other parent
                end
            else
                m_remaining(rem_temp_ordered(i_rem)) = m_remaining(rem_temp_ordered(i_rem))-m_rem_sort(i_rem); % subtract mass of remnant object from remaining mass of corresponding parent         
            end
        end    
        m_rem = m_rem_sort(idx_m_sort); %reorder masses of remnant objects
        rem1_temp = rem_temp_ordered(idx_m_sort)==1; %0s or 1s, 1 for remnant objects that are assigned to larger satellite
        d_rem_approx = zeros(numel(m_rem),1); %preallocate approximate diameters of remnant objects
        d_rem_approx(rem1_temp)  = (m_rem(rem1_temp) ./ p1_mass * p1_radius^3).^(1/3) * 2; % use same density (mass / r^3) as larger satellite to get approximate d_rem <<<<<<<<<<
        d_rem_approx(~rem1_temp) = (m_rem(~rem1_temp) ./ p2_mass * p2_radius^3).^(1/3) * 2; % use same density (mass / r^3) as smaller satellite to get approximate d_rem <<<<<<<<<<     
        Am_rem = func_Am(d_rem_approx, p1_objclass); %compute A/m using approximate diameters
        A_rem = m_rem.*Am_rem; %compute area from mass and A/m
        d_rem = d_rem_approx; %use approximate diameter based on parent object
        %remove remnant objects that are too small
        idx_too_small = find(d_rem<LB & m_rem<M/1000); %identify remnant objects that are smaller than lower bound Lc or less than 0.1% M
        m_rem(idx_too_small) = []; %remove remnant objects that are too small from mass vector
        d_rem(idx_too_small) = []; %remove remnant objects that are too small from diameter vector
        A_rem(idx_too_small) = []; %remove remnant objects that are too small from area vector
        Am_rem(idx_too_small) = []; %remove remnant objects that are too small from A/m vector
        rem1_temp(idx_too_small) = []; %remove remnant objects that are too small from vector of 0s or 1s, 1 if assigned to larger satellite, and 0 if assigned to smaller satellite
        idx_rem1 = find(rem1_temp);  %identify indices of remnant objects assigned to larger satellite
        idx_rem2 = find(~rem1_temp); %identify indices of remnant objects assigned to smaller satellite

    else    % Non-catastrophic collision fragments are deposited from 1 mm upward 
            % until the total mass, M = Mp * v_i^2  is achieved.  The final fragment 
            % is deposited in a single massive fragment reminiscent of a cratered target mass.
        % Assign remnant mass > d > A > Am > dv        
        largeidx = (d > p2_radius*2 | m > p2_mass) & d < p1_radius*2; %0s or 1s, 1 if larger than smaller satellite and shorter than larger satellite (except for remnant objects that have already been assigned to one of the satellites)
        m_assigned_large = max([0,sum(m(largeidx))]); %mass of fragments that need to belong to larger satellite
        
        if m_assigned_large>p1_mass %if assigned fragments to larger satellite already amount to more than the mass of the satellite
            % Sort by Lc, keep smallest objects until mass adds up <<<<<<<<<<<<<<
            idx_large = find(largeidx); %find indices of fragments assigned to larger satellite
            [~, dord1] = sort(m(idx_large));
            cumsum_m1 = cumsum(m(idx_large(dord1))); %cumulative sum of fragment masses in ascending order of diameter
            lastidx1 = find(cumsum_m1 < p1_mass, 1, 'last'); % last valid index of sorted list
            m(idx_large(dord1(lastidx1+1:end))) = [];   %Remove fragments exceeding mass constraint, that cannot belong to the smaller parent either
            d(idx_large(dord1(lastidx1+1:end))) = [];
            A(idx_large(dord1(lastidx1+1:end))) = [];
            Am(idx_large(dord1(lastidx1+1:end))) = [];
            largeidx(idx_large(dord1(lastidx1+1:end))) = [];
            m_assigned_large = cumsum_m1(lastidx1); %redefine mass assigned to larger satellite
            if isempty(m_assigned_large)  % NEW FIX
                m_assigned_large = 0; 
            end
        end        
        m_remaining_large = p1_mass-m_assigned_large; %remaining mass in larger satellite
        smallidx = d > p1_radius*2 & ~largeidx; %0s or 1s, 1 if larger than heavier satellite (this is a special case where sat2 is lighter but larger than sat1, and so fragments larger than sat1 are assigned to sat2)

        m_remSum = M - sum(m);                  % remaining mass
        if m_remSum>m_remaining_large %if remaining mass is larger than remaining mass of larger satellite, split mass between two parent satellites
            m_rem1 = m_remaining_large; %remnant object assigned to larger satellite
            d_rem_approx1 = (m_rem1 ./ p1_mass * p1_radius^3).^(1/3) * 2; % use same density to approximate diameter <<<<   
            m_rem2 = m_remSum-m_remaining_large; %remnant object assigned to smaller satellite
            d_rem_approx2 = (m_rem2 ./ p2_mass * p2_radius^3).^(1/3) * 2; % use same density to approximate diameter <<<<   
            d_rem_approx = [d_rem_approx1;d_rem_approx2]; %concatenate remnant objects assigned to each parent
            m_rem = [m_rem1;m_rem2]; %concatenate remnant objects assigned to each parent
            Am_rem = func_Am(d_rem_approx, p1_objclass); %compute A/m using approximate diameter
            A_rem = m_rem.*Am_rem; %compute area from mass and A/m
            d_rem = d_rem_approx; %use approximate diameter based on parent object
            idx_rem1 = 1; %index of remnant object assigned to larger satellite
            idx_rem2 = 2; %index of remnant object assigned to smaller satellite
        else %if remaining mass is smaller than remaining mass of larger satellite, assign remnant object to larger satellite
            m_rem = m_remSum;
            d_rem_approx = (m_rem ./ p1_mass * p1_radius^3).^(1/3) * 2; % use same density to approximate diameter <<<<   
            Am_rem = func_Am(d_rem_approx, p1_objclass); %compute A/m using approximate diameter
            A_rem = m_rem*Am_rem; %compute area from mass and A/m
            d_rem = d_rem_approx; %use approximate diameter based on parent object
            idx_rem1 = 1; %remnant object is always assigned to larger satellite
            idx_rem2 = []; %no remnant objects are assigned to smaller satellite
        end
        %remove remnant objects that are too small
        idx_too_small = find(d_rem<LB & m_rem<M/1000); %identify remnant objects that are smaller than lower bound Lc or less than 0.1% M
        m_rem(idx_too_small) = []; %remove remnant objects that are too small from mass vector
        d_rem(idx_too_small) = []; %remove remnant objects that are too small from diameter vector
        A_rem(idx_too_small) = []; %remove remnant objects that are too small from area vector
        Am_rem(idx_too_small) = []; %remove remnant objects that are too small from A/m vector
        if idx_too_small==1
            idx_rem1 = [];  %remove indices of remnant objects assigned to larger satellite
        elseif idx_too_small==2
            idx_rem2 = []; %remove indices of remnant objects assigned to smaller satellite
        end
    end

else
    % Sort by Lc, keep smallest objects until mass adds up <<<<<<<<<<<<<<
    [~, dord] = sort(m);
    cumsum_m = cumsum(m(dord)); %cumulative sum of fragment masses in ascending order of diameter
    lastidx = find(cumsum_m < M, 1, 'last'); % last valid index of sorted list
    cumsum_lower = max([0,cumsum_m(lastidx)]); %cumulative mass of ordered fragments that remain below M
    valididx = dord(1:lastidx); %indices of fragments fulfilling mass constraint    
    m = m(valididx);   %Select elements fulfilling mass constraint
    d = d(valididx);
    A = A(valididx);
    Am = Am(valididx); 
    largeidx = (d > p2_radius*2 | m > p2_mass) & d < p1_radius*2; %0s or 1s, 1 if larger than smaller satellite and shorter than larger satellite (except for remnant objects that have already been assigned to one of the satellites)
    smallidx = d > p1_radius*2 & ~largeidx; %0s or 1s, 1 if larger than heavier satellite (this is a special case where sat2 is lighter but larger than sat1, and so fragments larger than sat1 are assigned to sat2)
        
    %Check if there is mass remaining, and generate an additional fragment if needed
    m_rem = M-cumsum_lower; %remaining mass to accumulate a total mass of M
    if m_rem>M/1000 %if the remaining mass is larger than 0.1% M, add one more fragment (unless it is smaller than Lc)
        if m_rem>(p2_mass-max([0,sum(m(smallidx))])) %if mass of remaining fragment is larger than mass of smaller satellite that hasn't been filled yet, assign fragment to larger satellite
            rand_assign_frag = 1; %assign fragment to one of the satellites, 1 for larger satellite
        elseif m_rem>(p1_mass-max([0,sum(m(largeidx))])) %if mass of remaining fragment is larger than mass of larger satellite that hasn't been filled yet, assign fragment to smaller satellite
            rand_assign_frag = 2; %assign fragment to one of the satellites, 1 for larger satellite
        else %if remaining fragment could belong to either of the satellites
            rand_assign_frag = 1+round(rand(1)); %1 or 2, randomly assign fragment to one of the satellites, 1 for larger satellite, 2 for smaller satellite
        end
        switch rand_assign_frag %depending on which satellite the remaining fragment corresponds to
            case 1 %assign fragment to larger satellite
                d_rem_approx = (m_rem ./ p1_mass * p1_radius^3).^(1/3) * 2; % use same density (mass / r^3) as larger satellite to get approximate d_rem <<<<<<<<<<
                idx_rem1 = 1; idx_rem2 = [];
            case 2 %assign fragment to smaller satellite
                d_rem_approx = (m_rem ./ p2_mass * p2_radius^3).^(1/3) * 2; % use same density (mass / r^3) as smaller satellite to get approximate d_rem <<<<<<<<<<
                idx_rem1 = []; idx_rem2 = 1;
        end
        Am_rem = func_Am(d_rem_approx, p1_objclass); %compute A/m using approximate diameter
        A_rem = m_rem*Am_rem; %compute area from mass and A/m
        d_rem = d_rem_approx; %use approximate diameter based on parent object
        if d_rem<LB && m_rem<M/1000 %if remaining fragment is smaller than lower bound Lc, do not generate an additional fragment
            d_rem = []; A_rem = []; Am_rem = []; m_rem = []; idx_rem1 = []; idx_rem2 = [];  % no "remnants" exist
        end
    else %if there is no need for an additional fragment, do not generate one       
        d_rem = []; A_rem = []; Am_rem = []; m_rem = []; idx_rem1 = []; idx_rem2 = [];  % no "remnants" exist
    end
end
% figure
% scatter(d,m,'kx')
% hold on
% scatter(d_rem(idx_rem1),m_rem(idx_rem1),'rx')
% scatter(d_rem(idx_rem2),m_rem(idx_rem2),'bx')
% scatter(p1_radius*2,p1_mass,'ro')
% scatter(p2_radius*2,p2_mass,'bo')
% legend('General frags.','p1 frags.','p2 frags.','p1','p2','Location','northwest')
% % plot([d;d_rem],0.556945*[d;d_rem].^2.0047077./10.^(-0.925),'k-')
% set(gca,'XScale','log','YScale','log')
% title('m vs d')
% xlabel('d')
% ylabel('m')
% 
% figure
% scatter(d,Am,'kx')
% hold on
% scatter(d_rem(idx_rem1),Am_rem(idx_rem1),'rx')
% scatter(d_rem(idx_rem2),Am_rem(idx_rem2),'bx')
% scatter(p1_radius*2,pi*p1_radius^2/p1_mass,'ro')
% scatter(p1_radius*2,(2*p1_radius)^2/p1_mass,'rs')
% scatter(p2_radius*2,pi*p2_radius^2/p2_mass,'bo')
% scatter(p2_radius*2,(2*p2_radius)^2/p2_mass,'bs')
% legend('General frags.','p1 frags.','p2 frags.','p1 sphere','p1 cube','p2 sphere','p2 cube')
% set(gca,'XScale','log','YScale','log')
% title('A/m vs d')
% xlabel('d')
% ylabel('A/m')
% 
% figure
% scatter(d.^3,m,'kx')
% hold on
% scatter(d_rem(idx_rem1).^3,m_rem(idx_rem1),'rx')
% scatter(d_rem(idx_rem2).^3,m_rem(idx_rem2),'bx')
% scatter((p1_radius*2)^3,p1_mass,'ro')
% scatter((p2_radius*2)^3,p2_mass,'bo')
% legend('General frags.','p1 frags.','p2 frags.','p1','p2','Location','northwest')
% set(gca,'XScale','log','YScale','log')
% title('m vs d^3')
% xlabel('d^3')
% ylabel('m')

% Assign dv to random directions; create samples on unit sphere
dv = func_dv([Am; Am_rem],'col')/1000;       % km/s
u=rand(length(dv),1)*2-1;
theta=rand(length(dv),1)*2*pi;

v = sqrt(1 - u.^2);
p = [v.*cos(theta) v.*sin(theta) u];
try
    dv_vec = [p(:,1).*dv  p(:,2).*dv  p(:,3).*dv];
catch
    return
end

% create fragments
% NON-CATASTROPHIC: PUT IN PARENT OBJECT AS DEBRIS
if ~isCatastrophic
    % find remnant mass of parent object; add to *_rem object
    p1remnantmass = p1_mass + p2_mass - sum([m; m_rem]); % could be larger than p1_mass..!
    m_rem = [m_rem; p1remnantmass];
    d_rem = [d_rem; p1_radius*2];
    A_rem = [A_rem; pi*p1_radius^2]; % area not used for func_create_tlesv2_vec
    Am_rem = [Am_rem; pi*p1_radius^2 ./ p1remnantmass];   % AMR used only for Bstar in func_create_tlesv2_vec
    dv = [dv; 0];   % no change from parent object
    dv_vec = [dv_vec; [0,0,0]];  
end

if abs(sum([m; m_rem]) - p1_mass - p2_mass) > M*0.05  % FIXED TO USE SUM OF P1,P2 MASS, NOT M
    warning('Total sum of debris mass (%0.1f kg) differs from "mass" of original objects (%0.1f kg)', ...
        sum([m; m_rem]), p1_mass + p2_mass);
end

fragments = [[d; d_rem] [A; A_rem] [Am; Am_rem] [m; m_rem] dv dv_vec(:,1) dv_vec(:,2) dv_vec(:,3)];
% DEBUG: HISTOGRAM OF THESE ^
    % figure(10);clf;
    % subplot(511); histogram(d,100); xlabel('d (m)');
    % subplot(512); histogram(A,100); xlabel('A (m^2)');
    % subplot(513); histogram(Am,100); xlabel('Am (m^2/kg)');
    % subplot(514); histogram(m,100); xlabel('m (kg)');
    % subplot(515); histogram(dv * 1000,[0:10:1000]); xlabel('dv (m/s)');
    % subplot(511); title('fragment distributions', sprintf('Original mass: %0.1f kg, radius: %0.1f m',p1_mass, p1_radius));
%     pause;

% Attribute debris to p1 or p2

    % From the ESA's C++ code: 
    % First: All fragments larger than the smaller satellite (characteristic length) 
    % are assigned to the larger satellite (keeping track of cumulative mass).
    % Second: Assign fragments to the larger satellite until the cumulative mass 
    % equals the mass of the initial mass of the larger satellite. Fragments are 
    % randomly ordered, and so this assigns fragments of all sizes to the larger satellite.
    % Third: Assign remaining fragments to smaller satellite (still randomly 
    % ordered, and so these fragments are of any size too).
    
    % % index of debris too large to be from smaller object: assign to debris1
    % largeidx = find( fragments(:,1) > p2_radius*2 | fragments(:,4) > p2_mass );
    % largeidx = unique([largeidx; size(fragments,1)]); % add in last element to "large" col
    % smallidx = setdiff(1:size(fragments,1), largeidx);
    % smallidx = randperm(numel(smallidx));
    % % randomly assign to debris1 until m1 is filled
    % p1indx = find(sum(fragments(largeidx,4)) + cumsum(fragments(smallidx,4)) < p1_mass, 1, 'last');
    % largeidx = union(largeidx, smallidx(1:p1indx));
    % sum(fragments(largeidx,4))   % should be similar to p1_mass
    % smallidx = setdiff(1:size(fragments,1), largeidx) ;
    % sum(fragments(smallidx,4))   % should be similar to p2_mass

%Distribute fragments amongst parent 1 and parent 2
if isCatastrophic
    largeidx(end+1:end+numel(m_rem)) = 0; %concatenate placeholder for remnant objects (if any), which have already been assigned to one of the two satellites (if no remnant objects, nothing happens)
    largeidx(end-numel(m_rem)+idx_rem1) = 1; % add remnant objects assigned larger satellite
    smallidx(end+1:end+numel(m_rem)) = 0; %concatenate placeholder for remnant objects (if any), which have already been assigned to one of the two satellites (if no remnant objects, nothing happens)
    smallidx(end-numel(m_rem)+idx_rem2) = 1; % add remnant objects assigned smaller satellite
    if size(largeidx,2)>1 %make sure largeidx is a column vector
        largeidx = largeidx';
    end
    if size(smallidx,2)>1 %make sure smallidx is a column vector
        smallidx = smallidx';
    end
    assignedidx = largeidx | smallidx; %0s or 1s, 1 for fragments that have been assigned to one of the satellites so far
    idx_unassigned = find(~assignedidx); %indices of fragments that don't belong to any satellite yet
    % assign fragments to larger satellite until m1 is filled (fragments are randomly ordered
    % already): cumulative sum of fragments that have not been assigned to any satellite yet
    p1indx = max([0,find(sum(fragments(largeidx,4)) + cumsum(fragments(idx_unassigned,4)) <= p1_mass, 1, 'last')]);
    fragments1 = [fragments(largeidx,:); fragments(idx_unassigned(1:p1indx),:)];
    fragments2 = [fragments(smallidx,:); fragments(idx_unassigned(p1indx+1:end),:)];
else % non catastrophic: largest is with large object, rest is with smaller
    heaviestInd = fragments(:,4) == max(fragments(:,4)); % logical index, most massive
    lighterInd = ~(fragments(:,4) == max(fragments(:,4))); % logical index
    fragments1 = fragments(heaviestInd,:);
    fragments2 = fragments(lighterInd,:);
end

if any(fragments1(:,1) > p1_radius*2*1.00001)
    warning('Some fragments (%0.4f m) exceed the diameter of the larger object (%0.4f m)', ...
        max(fragments1(:,1)), p1_radius*2);
end

if any(fragments2(:,1) > p2_radius*2*1.00001)
    warning('Some fragments (%0.4f m) exceed the diameter of the smaller object (%0.4f m)', ...
        max(fragments2(:,1)), p2_radius*2);
end

%To match mass budget, some fragments smaller than LB may have been created
%>>> remove from list of fragments, after which, mass budget may not be
%fulfilled
fragments1(fragments1(:,1)<LB,:) = [];
fragments2(fragments2(:,1)<LB,:) = [];

% Create debris objects and add radius & mass directly
% mat_frag = [a,ecco,inclo,nodeo,argpo,mo,bstar,mass,radius,...
%             errors,controlled,a_desired,missionlife,constel,...
%             date_created,launch_date,r,v,objclass,ID_frag]; [num_sat x 24];
[debris1] = func_create_tlesv2_vec(ep, p1_r,p1_v,p1_objclass, fragments1, param.max_frag,param.mu,param.req,param.maxID); 
param.maxID = param.maxID+size(debris1,1);
[debris2] = func_create_tlesv2_vec(ep, p2_r,p2_v,p2_objclass, fragments2, param.max_frag,param.mu,param.req,param.maxID);

