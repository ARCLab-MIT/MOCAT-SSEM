function ceq =  equilib(x,data)

Dt = data.Dt;
N_shell = data.N_shell;

S_all = x(1:N_shell);
D_all = x(N_shell+1:2*N_shell);
N_all = x(2*N_shell+1:3*N_shell);
lambda_all = x(3*N_shell+1:4*N_shell); 

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

for k1=1:N_shell
    S_b(k1,1:N_shell-1) = S_all([1:k1-1, k1+1:end]);
    D_b(k1,1:N_shell-1) = D_all([1:k1-1, k1+1:end]);
    N_b(k1,1:N_shell-1) = N_all([1:k1-1, k1+1:end]);
end

n0_all_b = data.n0_all_b;
n1_all_b = data.n1_all_b;
n2_all_b = data.n2_all_b;
n3_all_b = data.n3_all_b;
n4_all_b = data.n4_all_b;
n5_all_b = data.n5_all_b;


yeq_S = zeros(1,N_shell);
yeq_D = yeq_S;
yeq_N = yeq_S;
for k = N_shell:-1:1

    S = S_all(k);
    D = D_all(k);
    N = N_all(k);
    lambda = lambda_all(k); 

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
    rvel_current_D = -rhok*beta(2)*sqrt(mu*(R0(k)))*(24*3600*365.25);
    rvel_current_N  =-rhok*beta(3)*sqrt(mu*(R0(k)))*(24*3600*365.25);

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

    if k == N_shell
        D_upper = 0; 
        N_upper = 0;
        % D_upper = D_all(k);  
        % N_upper = N_all(k);
    else
        D_upper = D_all(k+1); 
        N_upper = N_all(k+1); 
    end

    % Equality constraints: equilibrium equations
    yeq_S(k) = lambda + s1_all(k)*S + s2_all(k)*S*N + s3_all(k)*S*D + s4_all(k)*S^2;
    yeq_D(k) = d0_all(k)*S + d1_all(k)*D*S + d2_all(k)*S*N + d3_all(k)*D*N + d4_all(k)*D^2 + d5_all(k)*D_upper + d6_all(k)*D;
    yeq_N(k) = n0_all(k)*N*S + n1_all(k)*D*S + n2_all(k)*N*D + n3_all(k)*D^2 + n4_all(k)*S^2 + n5_all(k)*N^2 + n6_all(k)*N_upper + n7_all(k)*N;

    % modification to make N at one shell dependent on all the other ones
    for k1=N_shell-1:-1:1
        yeq_N_b(k,k1) = n0_all_b(k,k1)*N_b(k,k1)*S_b(k,k1) + ...
                        n1_all_b(k,k1)*D_b(k,k1)*S_b(k,k1) + ...
                        n2_all_b(k,k1)*N_b(k,k1)*D_b(k,k1) + ...
                        n3_all_b(k,k1)*D_b(k,k1)^2 + ...
                        n4_all_b(k,k1)*S_b(k,k1)^2 + ...
                        n5_all_b(k,k1)*N_b(k,k1)^2;
    end
    yeq_N(k) = yeq_N(k) + sum(yeq_N_b(k,:),2);

end

ceq = [yeq_S , yeq_D , yeq_N]';

end