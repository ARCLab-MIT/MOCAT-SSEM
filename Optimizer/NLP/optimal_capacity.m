close all
clear 
clc

%==========================================================================
%==========================================================================
% Description:

% A constrained nonlinear programming (NLP) optimization problem is 
% formulated to compute the LEO optimal orbital carrying capacity solutions, 
% considering stable equilibrium points and sustainability constraints. 
% A three species MOCAT-SSEM is considered, with active satellites (S), 
% derelicts (D), and debris (N).

% Author: Dr. Giovanni Lavezzi 
% v1: 08/2024
% v2: 03/2025
%==========================================================================
%==========================================================================
% Data

sel_plots = 1; % plot selection, 1: yes, 0: no

N_shell = 24;
h_min = 200; % initial altitude [km]
h_max = 1400; % final altitude [km] 2000
R0 = linspace(h_min,h_max,N_shell+1);
R02 = linspace(h_min,h_max,N_shell+1);
deltaH = R02(2)-R02(1); % thickness of the shell [km]
re = 6378.1; % [km]
years_s = 365*24*3600;

mu = 3.986004418e14;%meters^3/s^2
R1 = R0;
R0 = (re+R0)*1000;
V = 4/3*pi*(diff(R0.^3)); % volume of the shells [m^3]

Dt = 8; % operation lifetime: after those years, active satellites become derelict [years] 8
v = 10*1000*(24*3600*365.25);% impact velocity [m/year]
Dhl = deltaH*1000;
Dhu = -deltaH*1000;
Cd = 2.2;
failure_rate_L = 0; % lower boundary
failure_rate_U = 0.1; %3.5e-4; % upper boundary 

failure_rate_L = failure_rate_L/100;
failure_rate_U = failure_rate_U/100; % [% = fail_rate / 100] 

P = 0.9;                
alpha = 0.1; % fraction of collisions that satellites fail to avoid 0.2      
alpha_active = 0.01; % fraction of collisions occuring between active satellites 0.01             
delta = 0;

%==========================================================================
% physical characteristics

mass      = [200.0  200.0   0.5];
diameter  = [2.5    2.5    0.25];
Area      = [5.0    5.0    0.05];

n_obj = length(mass);
radii = diameter/2;
area_mass = Area./mass; % area/mass ratio  [# shells, # object types]
beta = Cd*area_mass; % ballistic coefficient [# shells, # object types]

%==========================================================================
% NASA standard breakup models implementation for collisions
v_imp = 10; % impact velocity [km/s]
LC = 0.1; % minimum size of fragments [m]
n_f_catastrophic = @(M1,M2) 0.1*LC^(-1.71)*(M1+M2)^(0.75); % number of fragments generated during a catastrophic collision (NASA standard break-up model). M is the sum of the mass of the objects colliding in kg
n_f_damaging = @(M1,M2) 0.1*LC^(-1.71)*(min(M1,M2)*v_imp^2)^(0.75); % number of fragments generated during a non-catastrophic collision (improved NASA standard break-up model: takes into account the kinetic energy). M is the mass of the less massive object colliding in kg
ind_obj_catastrophic = inf; % all collisions are catastrophic
sigma = zeros(n_obj,n_obj);
for q = 1:length(radii)
    mass_q = mass(q);
    for qq = 1:length(radii)
        sigma(q,qq) = (radii(q)+radii(qq))^2;
        mass_qq = mass(qq);
    end
end
f_corrective = 1; 
phi = f_corrective*pi*sigma*v; 

% number of fragments generated in each collisions among the species for each altitude bin
[n1,~,~,hc1] = EVOLVEbinsDV(mass(1),mass(1),radii(1),radii(1),v_imp,[],[mass(3),inf],[],LC,0,R02);  % S-S and S-D and D-D
[n2,~,~,hc2] = EVOLVEbinsDV(mass(1),mass(3),radii(1),radii(3),v_imp,[],[mass(3),inf],[],LC,0,R02);  % S-N and D-N
[n3,~,~,hc3] = EVOLVEbinsDV(mass(3),mass(3),radii(3),radii(3),v_imp,[],[mass(3),inf],[],LC,0,R02);  % N-N

K0(1:2,1:2) = n1;
K0(1:2,3) = n2;
K0(3,1:2) = n2;
K0(3,3) = n3;

data.hc1 = hc1;
data.hc2 = hc2;
data.hc3 = hc3;

%==========================================================================
%% Optimization

fprintf("Optimization:\n")
data.failure_rate_L = failure_rate_L;
data.failure_rate_U = failure_rate_U;
data.Dt = Dt;
data.N_shell = N_shell;

data.mu = mu;
data.phi = phi;
data.V = V;
data.K0 = K0;
data.R0 = R0;
data.R02 = R02;
data.Dhu = Dhu;
data.Dhl = Dhl;

data.alpha = alpha;
data.beta = beta;
data.delta = delta;
data.P = P;
data.alpha_active = alpha_active;

x0_S = ones(1,N_shell)*1e3 *0;
x0_D = ones(1,N_shell)*1e3 *0;
x0_N = ones(1,N_shell)*1e3 *0;

x0_lambda = ones(1,N_shell)*1e3 *0;

x0_S = x0_S + 0;
x0_D = x0_D + 0; 
x0_N = x0_N + 0;
x0_lambda = x0_lambda + 0;
x0 = [x0_S, x0_D, x0_N, x0_lambda] + 0;
options = optimoptions( 'fmincon' , 'Display' , 'off' ,...
        'Algorithm', 'sqp', ...
        'MaxIterations', 5e5*1e0, ... 
        'MaxFunctionEvaluations', 5e6*1e0, ... 
        'ConstraintTolerance', 1e-11, ...
        'FunctionTolerance', 1e-12, ... 
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-18, ... 
        'UseParallel', false );
A = []; b = []; Aeq = []; beq = []; 

lb_S = 0.45 * ones(size(x0_S));
lb_D = 0.0005 * ones(size(x0_D));
lb_N = 0.0005 * ones(size(x0_N));
lb_l = 0.0005 * ones(size(x0_lambda)); 

max_shell = 14 + 1;
lb_S(max_shell:end) = 0;
lb_D(max_shell:end) = 0;
lb_N(max_shell:end) = 0;
lb_l(max_shell:end) = 0;

lb = [lb_S, ...
      lb_D, ...
      lb_N, ...
      lb_l];
ub = [];
    
[N_shell,...
s1_all,s2_all,s3_all,s4_all,...
d0_all,d1_all,d2_all,d3_all,d4_all,d5_all,d6_all, ...
n0_all,n1_all,n2_all,n3_all,n4_all,n5_all,n6_all,n7_all, ...
n0_all_b, n1_all_b, n2_all_b, n3_all_b, n4_all_b, n5_all_b, ...
Dt,failure_rate_L,failure_rate_U] = data_for_solver(data);

data.n0_all_b = n0_all_b; 
data.n1_all_b = n1_all_b;
data.n2_all_b = n2_all_b; 
data.n3_all_b = n3_all_b; 
data.n4_all_b = n4_all_b; 
data.n5_all_b = n5_all_b;

nonlcon = @(x) PrcAll(x,N_shell,...
s1_all,s2_all,s3_all,s4_all,...
d0_all,d1_all,d2_all,d3_all,d4_all,d5_all,d6_all, ...
n0_all,n1_all,n2_all,n3_all,n4_all,n5_all,n6_all,n7_all, ...
n0_all_b, n1_all_b, n2_all_b, n3_all_b, n4_all_b, n5_all_b, ...
Dt,failure_rate_L,failure_rate_U);
objective = @(x) PrfAll(x,data);

%======================================================================
cpu_i = tic;
[xopt,fval,exitflag,output] = fmincon( objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, options );
fval
cpu_f = toc(cpu_i);
if exitflag < 1 || exitflag == 2
    output.message
end

%======================================================================
% Define output variables
X_out(1,:) = xopt(1:N_shell);
X_out(2,:) = xopt(N_shell+1:2*N_shell);
X_out(3,:) = xopt(2*N_shell+1:3*N_shell);
lam1(1,:) = xopt(3*N_shell+1:4*N_shell); 
for i1 = 1:N_shell
    if X_out(1,i1) < 1e-2
        X_out(1,i1) = 0;
    end
    if X_out(2,i1) < 1e-2
        X_out(2,i1) = 0;
    end
    if X_out(3,i1) < 1e-2
        X_out(3,i1) = 0;
    end
    if lam1(1,i1) < 1e-2
        lam1(1,i1) = 0;
    end
end

fprintf("cpu time: %f s\n",cpu_f)
fprintf("N_shell: %d \n",N_shell)
fprintf("Altitude range: %d-%d km\n",h_min,h_max)
fprintf("Constraint Tolerance: %e\n",options.ConstraintTolerance)
fprintf("Function Tolerance: %e\n",options.FunctionTolerance)

ceq = equilib(xopt,data);
S_eq = ceq(1:N_shell);
D_eq = ceq(N_shell+1:2*N_shell);
N_eq = ceq(2*N_shell+1:3*N_shell);


%% Verification 

S = sym('S',[N_shell,1]);
D = sym('D',[N_shell,1]);
N = sym('N',[N_shell,1]);
lam = sym('lam',[N_shell,1]);

for k1=1:N_shell
    S_b(k1,1:N_shell-1) = S([1:k1-1, k1+1:end]);
    D_b(k1,1:N_shell-1) = D([1:k1-1, k1+1:end]);
    N_b(k1,1:N_shell-1) = N([1:k1-1, k1+1:end]);

    K0_SS_b(k1,1:N_shell-1) = [hc1(N_shell-k1+1:N_shell-1), flip(hc1(end-k1+1:-1:N_shell+1))];
    K0_SD_b(k1,1:N_shell-1) = [hc1(N_shell-k1+1:N_shell-1), flip(hc1(end-k1+1:-1:N_shell+1))];
    K0_SN_b(k1,1:N_shell-1) = [hc2(N_shell-k1+1:N_shell-1), flip(hc2(end-k1+1:-1:N_shell+1))];
    K0_DD_b(k1,1:N_shell-1) = [hc1(N_shell-k1+1:N_shell-1), flip(hc1(end-k1+1:-1:N_shell+1))];
    K0_DN_b(k1,1:N_shell-1) = [hc2(N_shell-k1+1:N_shell-1), flip(hc2(end-k1+1:-1:N_shell+1))];
    K0_NN_b(k1,1:N_shell-1) = [hc3(N_shell-k1+1:N_shell-1), flip(hc3(end-k1+1:-1:N_shell+1))];
end

for k=1:N_shell 
    
    if k<N_shell
        n0=(N(k+1));
        n_upper=D(k+1);
        rho_k1 = densityexp(R02(k+1)); % no solar flux
        rvel_upper_D=-rho_k1*beta(2)*sqrt(mu*(R0(k+1)))*(24*3600*365.25);% drag flux
        rvel_upper_N=-rho_k1*beta(3)*sqrt(mu*(R0(k+1)))*(24*3600*365.25);%drag flux
        
    else
        n_upper = 0;
        n0 = 0;
        rho_k1 = densityexp((R02(k+1)+(R02(k+1)-R02(k)))); % no solar flux
        rvel_upper_D=-rho_k1*beta(2)*sqrt(mu*(R0(k+1)+R0(k+1)-R0(k)))*(24*3600*365.25);% drag flux
        rvel_upper_N=-rho_k1*beta(3)*sqrt(mu*(R0(k+1)+R0(k+1)-R0(k)))*(24*3600*365.25);%drag flux
    end
    rhok = densityexp(R02(k)); % no solar flux
    rvel_current_D=-rhok*beta(2)*sqrt(mu*(R0(k)))*(24*3600*365.25);
    rvel_current_N=-rhok*beta(3)*sqrt(mu*(R0(k)))*(24*3600*365.25);

    phi_SS = phi(1,1)/V(k);
    phi_SD = phi(1,2)/V(k);
    phi_SN = phi(1,3)/V(k);
    phi_DD = phi(2,2)/V(k);
    phi_DN = phi(2,3)/V(k);
    phi_NN = phi(3,3)/V(k); 
    K0_SS = K0(1,1);
    K0_SD = K0(1,2);
    K0_SN = K0(1,3);
    K0_DD = K0(2,2);
    K0_DN = K0(2,3);
    K0_NN = K0(3,3);

    Eq1(k)=lam(k)-S(k)/Dt-(delta+alpha)*(phi_SN*N(k)*S(k)+phi_SD*D(k)*S(k))...
        -alpha_active*phi_SS*S(k)^2;

    Eq2(k)=(1-P)*S(k)/Dt+...
        delta*phi_SD*D(k)*S(k)+delta*phi_SN*N(k)*S(k)-(phi_DN*N(k)*D(k)+phi_DD*D(k)*D(k))...
        +n_upper*rvel_upper_D/Dhu+D(k)*rvel_current_D/Dhl;%flux

    Eq3(k)=(K0_SN*phi_SN*N(k)*(alpha*S(k))+K0_SD*phi_SD*D(k)*(alpha*S(k)))+(K0_DN*phi_DN*N(k)*D(k)+...
        K0_DD*phi_DD*D(k)*D(k))+alpha_active*K0_SS*phi_SS*S(k)^2+phi_NN*K0_NN*N(k)^2 ...
            +n0*rvel_upper_N/Dhu+N(k)*rvel_current_N/Dhl; 

    % modification to make N at one shell dependent on all the other ones
    for k1=1:N_shell-1
        Eq3b(k,k1) = K0_SN_b(k,k1)*phi_SN*N_b(k,k1)*(alpha*S_b(k,k1)) + ...
                   K0_SD_b(k,k1)*phi_SD*D_b(k,k1)*(alpha*S_b(k,k1)) + ...
                   K0_DN_b(k,k1)*phi_DN*N_b(k,k1)*D_b(k,k1) + ...
                   K0_DD_b(k,k1)*phi_DD*D_b(k,k1)*D_b(k,k1) + ...
                   alpha_active*K0_SS_b(k,k1)*phi_SS*S_b(k,k1)^2 + ...
                   phi_NN*K0_NN_b(k,k1)*N_b(k,k1)^2;
    end
    Eq3(k) = Eq3(k) + sum(Eq3b(k,:),2);

end

var=[S;D;N];
f3=[Eq1.';Eq2.';Eq3.'];

fun3 = matlabFunction(f3,'Vars',{var,lam});
func_temp = @(x) (fun3(x,lam1'));
func = @(t,x) func_temp(x);

tf_ss = 1000;
tspan1 = linspace(0,tf_ss,1000);
options_ode = odeset('reltol', 1e-8,'abstol', 1e-8);
[t_prop,x_prop] = ode15s(func, tspan1, xopt(1:3*N_shell), options_ode);

S_prop = x_prop(:,1:N_shell);
D_prop = x_prop(:,N_shell+1:2*N_shell);
N_prop = x_prop(:,2*N_shell+1:3*N_shell);
N_tot = S_prop + D_prop + N_prop;
N_tot_sum = sum(N_tot,2);
S_sum = sum(S_prop,2);
D_sum = sum(D_prop,2);
N_sum = sum(N_prop,2);

%==========================================================================
%% Plots

if sel_plots 
    
    colors = {[0, 0, 0], [0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]};
    sel_LineWidth = 2;
    sel_MarkerWidth = 10;
    sel_LineWidthAxis = 1;
    sel_FontSize = 14;
    
    figure('Color', 'w');
    hold on; grid on;
    plot(t_prop,N_tot_sum(end)-N_tot_sum,'Color',colors{1},'LineWidth',sel_LineWidth)
    plot(t_prop,S_sum(end)-S_sum,'Color',colors{2},'LineWidth',sel_LineWidth);
    plot(t_prop,D_sum(end)-D_sum,'Color',colors{3},'LineWidth',sel_LineWidth);
    plot(t_prop,N_sum(end)-N_sum,'Color',colors{4},'LineWidth',sel_LineWidth);
    xlabel("Years");ylabel('Count')
    legend("Total","S","D","N","Location","best")
    set(gca, 'FontSize', sel_FontSize, 'linewidth', sel_LineWidthAxis)
    
    figure('Color', 'w');
    hold on; grid on;
    plot(R02(2:end)- 25,100*failure_rate_U*ones(size(R02(2:end))),'--','Color',colors{1},'LineWidth',sel_LineWidth);
    plot(R02(2:end) - 25,100*(Dt*lam1-X_out(1,:))./(Dt*lam1),'-o','Color',colors{6},'LineWidth',sel_LineWidth,'MarkerSize',sel_MarkerWidth);
    xlabel('Altitude (km)');
    ylabel('Failure rate (%)');
    legend("\chi_{max}","Location","best")
    set(gca, 'FontSize', sel_FontSize, 'linewidth', sel_LineWidthAxis)
    
    figure('Color', 'w');
    semilogy(R02(2:end) - 25,X_out(1,:),'-o','Color',colors{2},'LineWidth',sel_LineWidth,'MarkerSize',sel_MarkerWidth);
    hold on; grid on;
    semilogy(R02(2:end) - 25,X_out(2,:),'-o','Color',colors{3},'LineWidth',sel_LineWidth,'MarkerSize',sel_MarkerWidth);
    semilogy(R02(2:end) - 25,X_out(3,:),'-o','Color',colors{4},'LineWidth',sel_LineWidth,'MarkerSize',sel_MarkerWidth);
    semilogy(R02(2:end) - 25,lam1(1,:),'-v','Color',colors{5},'LineWidth',sel_LineWidth,'MarkerSize',sel_MarkerWidth);
    xlabel('Altitude (km)');ylabel('Count');
    legend('S','D','N','\lambda',"Location","best")
    set(gca, 'FontSize', sel_FontSize, 'linewidth', sel_LineWidthAxis)
    
    figure('Color', 'w');
    hold on; grid on;
    plot(R02(2:end) - 25,S_eq,'Color',colors{2},'LineWidth',sel_LineWidth);
    plot(R02(2:end) - 25,D_eq,'Color',colors{3},'LineWidth',sel_LineWidth);
    plot(R02(2:end) - 25,N_eq,'Color',colors{4},'LineWidth',sel_LineWidth);
    xlabel("Altitude (km)");ylabel('Count')
    legend("S","D","N","Location","best")
    set(gca, 'FontSize', sel_FontSize, 'linewidth', sel_LineWidthAxis)

    % Stability plot
    xopt_mod = xopt(1:3*N_shell);
    lam1_mod(1,:) = xopt(3*N_shell+1:4*N_shell); 
    
    tf_ss = 1000;
    tspan1 = linspace(0,tf_ss,1000);
    options_ode = odeset('reltol', 1e-8,'abstol', 1e-8);

    factor_lambda = [];
    factor_init = 1; % [%]
    factor_end = 200; % [%]
    factor_step = 5;

    factor = (factor_init:factor_step:factor_end)/100; 
    factor_N = (factor_init:factor_step:factor_end)/100; 
    
    N_tot_sum_mod = []; N_tot_mod = []; 
    t_prop_mod = []; t_prop_mod_all = [];
    x_prop_mod_all = [];
    S_sum_mod = []; D_sum_mod = []; N_sum_mod = [];
    S_prop_mod = []; D_prop_mod = []; N_prop_mod = [];

    for ii1 = 1:length(factor)
    
        if length(factor_lambda) > 1
            func_temp2 = @(x) (fun3(x,lam1_mod' * factor_lambda(ii1)));
        elseif length(factor_lambda) == 1
            func_temp2 = @(x) (fun3(x,lam1_mod' * factor_lambda));
        else
            func_temp2 = @(x) (fun3(x,lam1_mod'));
        end
        func2 = @(t,x) func_temp2(x);
        [t_prop_mod_temp,x_prop_mod] = ode15s(func2, tspan1, [xopt_mod(:,1:N_shell) * factor(ii1) , xopt_mod(:,N_shell+1:2*N_shell) * factor(ii1), xopt_mod(:,2*N_shell+1:3*N_shell) * factor_N(ii1)], options_ode);
        
        t_prop_mod_all{ii1} = t_prop_mod_temp;
        x_prop_mod_all{ii1} = x_prop_mod;
    
        S_prop_mod{ii1} = x_prop_mod(:,1:N_shell);
        D_prop_mod{ii1} = x_prop_mod(:,N_shell+1:2*N_shell);
        N_prop_mod{ii1} = x_prop_mod(:,2*N_shell+1:3*N_shell);
        N_tot_mod{ii1} = S_prop_mod{ii1} + D_prop_mod{ii1} + N_prop_mod{ii1};
        N_tot_sum_mod(1:length(sum(N_tot_mod{ii1},2)),ii1) = sum(N_tot_mod{ii1},2);
        S_sum_mod(1:length(sum(S_prop_mod{ii1},2)),ii1) = sum(S_prop_mod{ii1},2);
        D_sum_mod(1:length(sum(D_prop_mod{ii1},2)),ii1) = sum(D_prop_mod{ii1},2);
        N_sum_mod(1:length(sum(N_prop_mod{ii1},2)),ii1) = sum(N_prop_mod{ii1},2);
        t_prop_mod(1:length(sum(t_prop_mod_temp,2)),ii1) = t_prop_mod_temp;
    
    end
    figure('Color', 'w');
    h11 = semilogy(t_prop,S_sum,'-','Color',colors{2},'LineWidth',sel_LineWidth,'MarkerSize',sel_MarkerWidth);
    hold on; grid on;
    h12 = semilogy(t_prop,D_sum,'-','Color',colors{3},'LineWidth',sel_LineWidth,'MarkerSize',sel_MarkerWidth);
    h13 = semilogy(t_prop,N_sum,'-','Color',colors{4},'LineWidth',sel_LineWidth,'MarkerSize',sel_MarkerWidth);
    h21 = semilogy(t_prop_mod,S_sum_mod,':','Color',colors{2},'LineWidth',sel_LineWidth,'MarkerSize',sel_MarkerWidth);
    h22 = semilogy(t_prop_mod,D_sum_mod,':','Color',colors{3},'LineWidth',sel_LineWidth,'MarkerSize',sel_MarkerWidth);
    h23 = semilogy(t_prop_mod,N_sum_mod,':','Color',colors{4},'LineWidth',sel_LineWidth,'MarkerSize',sel_MarkerWidth);
    xlabel("Years");ylabel('Count')
    xlim([0 1000]);
    legend([h11(1) h12(1) h13(1) h21(1) h22(1) h23(1)],{"S_{equil}","D_{equil}","N_{equil}","S","D","N"},"Location","best")
    set(gca, 'FontSize', sel_FontSize, 'linewidth', sel_LineWidthAxis)

end

%==========================================================================