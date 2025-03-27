function [c,ceq,dc,dceq] = PrcAll(x, N_shell,...
    s1_all,s2_all,s3_all,s4_all,...
    d0_all,d1_all,d2_all,d3_all,d4_all,d5_all,d6_all, ...
    n0_all,n1_all,n2_all,n3_all,n4_all,n5_all,n6_all,n7_all, ...
    n0_all_b, n1_all_b, n2_all_b, n3_all_b, n4_all_b, n5_all_b, ...
    Dt,failure_rate_L,failure_rate_U)

S_all = x(1:N_shell);
D_all = x(N_shell+1:2*N_shell);
N_all = x(2*N_shell+1:3*N_shell);
lambda_all = x(3*N_shell+1:4*N_shell); 

yeq_S = zeros(1,N_shell);
yeq_D = yeq_S;
yeq_N = yeq_S;

for k1=1:N_shell
    S_b(k1,1:N_shell-1) = S_all([1:k1-1, k1+1:end]);
    D_b(k1,1:N_shell-1) = D_all([1:k1-1, k1+1:end]);
    N_b(k1,1:N_shell-1) = N_all([1:k1-1, k1+1:end]);
end

for k=1:N_shell

    S = S_all(k);
    D = D_all(k);
    N = N_all(k);
    lambda = lambda_all(k); 

    if k == N_shell
        D_upper = 0; 
        N_upper = 0;
    else
        D_upper = D_all(k+1); 
        N_upper = N_all(k+1); 
    end

    % Equality constraints: equilibrium equations
    yeq_S(k) = lambda + s1_all(k)*S + s2_all(k)*S*N + s3_all(k)*S*D + s4_all(k)*S^2;
    yeq_D(k) = d0_all(k)*S + d1_all(k)*D*S + d2_all(k)*S*N + d3_all(k)*D*N + d4_all(k)*D^2 + d5_all(k)*D_upper + d6_all(k)*D;
    yeq_N(k) = n0_all(k)*N*S + n1_all(k)*D*S + n2_all(k)*N*D + n3_all(k)*D^2 + n4_all(k)*S^2 + n5_all(k)*N^2 + n6_all(k)*N_upper + n7_all(k)*N;

    % modification to make N at one shell dependent on all the other ones
    for k1=1:N_shell-1
        yeq_N_b(k,k1) = n0_all_b(k,k1)*N_b(k,k1)*S_b(k,k1) + ...
                        n1_all_b(k,k1)*D_b(k,k1)*S_b(k,k1) + ...
                        n2_all_b(k,k1)*N_b(k,k1)*D_b(k,k1) + ...
                        n3_all_b(k,k1)*D_b(k,k1)^2 + ...
                        n4_all_b(k,k1)*S_b(k,k1)^2 + ...
                        n5_all_b(k,k1)*N_b(k,k1)^2;
    end
    yeq_N(k) = yeq_N(k) + sum(yeq_N_b(k,:),2);

    % Inequality constraints - failure rate
    % y_fail_l(k) = -(lambda * Dt * (1 - failure_rate_L) - S);
    % y_fail_u(k) = lambda * Dt * (1 - failure_rate_U) - S; 

    % Inequality constraints - max debris and derelicts per shell
    y_fail_l(k) = N - 1500;
    y_fail_u(k) = D - 300; 

end

c = [y_fail_u , y_fail_l].'; 
ceq = [yeq_S , yeq_D , yeq_N]';

dc = [];
dceq = [];

end