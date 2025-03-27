function [N_shell,...
s1_all,s2_all,s3_all,s4_all,...
d0_all,d1_all,d2_all,d3_all,d4_all,d5_all,d6_all, ...
n0_all,n1_all,n2_all,n3_all,n4_all,n5_all,n6_all,n7_all, ...
n0_all_b, n1_all_b, n2_all_b, n3_all_b, n4_all_b, n5_all_b, ...
Dt,failure_rate_L,failure_rate_U] = data_for_solver(data)

failure_rate_L = data.failure_rate_L;
failure_rate_U = data.failure_rate_U;
Dt = data.Dt;
N_shell = data.N_shell;

mu = data.mu;
phi = data.phi;
V = data.V;
K0 = data.K0;
R0 = data.R0;
R02 = data.R02;
Dhu = data.Dhu;
Dhl = data.Dhl;

alpha = data.alpha;
beta = data.beta;
delta = data.delta;
P = data.P;
alpha_active = data.alpha_active;

hc1 = data.hc1;
hc2 = data.hc2;
hc3 = data.hc3;

for k1=1:N_shell
    K0_SS_b(k1,1:N_shell-1) = [hc1(N_shell-k1+1:N_shell-1), flip(hc1(end-k1+1:-1:N_shell+1))];
    K0_SD_b(k1,1:N_shell-1) = [hc1(N_shell-k1+1:N_shell-1), flip(hc1(end-k1+1:-1:N_shell+1))];
    K0_SN_b(k1,1:N_shell-1) = [hc2(N_shell-k1+1:N_shell-1), flip(hc2(end-k1+1:-1:N_shell+1))];
    K0_DD_b(k1,1:N_shell-1) = [hc1(N_shell-k1+1:N_shell-1), flip(hc1(end-k1+1:-1:N_shell+1))];
    K0_DN_b(k1,1:N_shell-1) = [hc2(N_shell-k1+1:N_shell-1), flip(hc2(end-k1+1:-1:N_shell+1))];
    K0_NN_b(k1,1:N_shell-1) = [hc3(N_shell-k1+1:N_shell-1), flip(hc3(end-k1+1:-1:N_shell+1))];
end

for k=1:N_shell
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

    rhok = densityexp(R02(k)); % no solar flux
    rvel_current_D =-rhok*beta(2)*sqrt(mu*(R0(k)))*(24*3600*365.25);
    rvel_current_N =-rhok*beta(3)*sqrt(mu*(R0(k)))*(24*3600*365.25);

    if k<N_shell
        rho_k1 = densityexp(R02(k+1)); % no solar flux
        rvel_upper_D=-rho_k1*beta(2)*sqrt(mu*(R0(k+1)))*(24*3600*365.25);% drag flux
        rvel_upper_N=-rho_k1*beta(3)*sqrt(mu*(R0(k+1)))*(24*3600*365.25);%drag flux
        
    else
        rho_k1 = densityexp((R02(k+1)+(R02(k+1)-R02(k)))); % no solar flux
        rvel_upper_D=-rho_k1*beta(2)*sqrt(mu*(R0(k+1)+R0(k+1)-R0(k)))*(24*3600*365.25);% drag flux
        rvel_upper_N=-rho_k1*beta(3)*sqrt(mu*(R0(k+1)+R0(k+1)-R0(k)))*(24*3600*365.25);%drag flux
    end

    s1_all(k) = -1/Dt;
    s2_all(k) = -(delta+alpha)*phi_SN;
    s3_all(k) = -(delta+alpha)*phi_SD;
    s4_all(k) = -alpha_active*phi_SS;

    d0_all(k) = (1-P)/Dt; % S(k)
    d1_all(k) = delta*phi_SD; % D(k)*S(k)
    d2_all(k) = delta*phi_SN; % (N(k))*S(k)
    d3_all(k) = -phi_DN; % N(k)*D(k)
    d4_all(k) = -phi_DD; % (D(k))*D(k)
    d5_all(k) = rvel_upper_D/Dhu; % D(k+1)
    d6_all(k) = rvel_current_D/Dhl; % D(k)

    n0_all(k) = K0_SN*phi_SN*alpha; % (N(k))*S(k)
    n1_all(k) = K0_SD*phi_SD*alpha; % D(k)*S(k)
    n2_all(k) = K0_DN*phi_DN; % (N(k))*D(k)
    n3_all(k) = K0_DD*phi_DD; % (D(k))*D(k)
    n4_all(k) = alpha_active*K0_SS*phi_SS; % S(k)^2
    n5_all(k) = phi_NN*K0_NN; % (N(k))^2
    n6_all(k) = rvel_upper_N/Dhu; % N(k+1)
    n7_all(k) = rvel_current_N/Dhl; % N(k)  

    % for k1=N_shell-1:-1:1
    for k1=1:N_shell-1
        n0_all_b(k,k1) = K0_SN_b(k,k1)*phi_SN*alpha; % (N(k))*S(k)
        n1_all_b(k,k1) = K0_SD_b(k,k1)*phi_SD*alpha; % D(k)*S(k)
        n2_all_b(k,k1) = K0_DN_b(k,k1)*phi_DN; % (N(k))*D(k)
        n3_all_b(k,k1) = K0_DD_b(k,k1)*phi_DD; % (D(k))*D(k)
        n4_all_b(k,k1) = alpha_active*K0_SS_b(k,k1)*phi_SS; % S(k)^2
        n5_all_b(k,k1) = phi_NN*K0_NN_b(k,k1); % (N(k))^2    
    end

end

end

