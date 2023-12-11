% Verification with hardcoded version
% this is for null launch rate, and using same x0 from workbook in the mocat3_eq_check.m
% try for different launch rate (use same x0 and launch rate in workbook and hardcoded version)
% try with static exponential model
function verify_equations(my_sim)
    
    if length(my_sim.species_list) == 3
    % S
    disp('S')
    vpa(my_sim.equations(:,1),4)
    
    % hardcoded from mocat3_eq_check
    disp('lam1 - 0.2*S1 - 9.067e-8*N1*S1 - 2.831e-10*S1^2 - 2.887e-7*D1*S1')
    disp('lam2 - 0.2*S2 - 8.697e-8*N2*S2 - 2.715e-10*S2^2 - 2.769e-7*D2*S2')
    disp('lam3 - 0.2*S3 - 8.349e-8*N3*S3 - 2.606e-10*S3^2 - 2.659e-7*D3*S3')
    disp('lam4 - 0.2*S4 - 8.022e-8*N4*S4 - 2.504e-10*S4^2 - 2.554e-7*D4*S4')
    disp('lam5 - 0.2*S5 - 7.713e-8*N5*S5 - 2.408e-10*S5^2 - 2.456e-7*D5*S5')
    
    
    % D
    disp('D')
    vpa(my_sim.equations(:,2),4)
    
    % hardcoded from mocat3_eq_check
    disp('2.298*D2 - 55.29*D1 + 0.01*S1 + 8.89e-8*N1*S1 - 2.831e-8*D1^2 - 8.89e-9*D1*N1 + 2.831e-7*D1*S1')
    disp('0.1959*D3 - 2.298*D2 + 0.01*S2 + 8.527e-8*N2*S2 - 2.715e-8*D2^2 - 8.527e-9*D2*N2 + 2.715e-7*D2*S2')
    disp('0.02251*D4 - 0.1959*D3 + 0.01*S3 + 8.186e-8*N3*S3 - 2.606e-8*D3^2 - 8.186e-9*D3*N3 + 2.606e-7*D3*S3')
    disp('0.003794*D5 - 0.02251*D4 + 0.01*S4 + 7.864e-8*N4*S4 - 2.504e-8*D4^2 - 7.864e-9*D4*N4 + 2.504e-7*D4*S4')
    disp('0.01*S5 - 0.003794*D5 + 7.562e-8*N5*S5 - 2.408e-8*D5^2 - 7.562e-9*D5*N5 + 2.408e-7*D5*S5')
    
    
    % N
    disp('N')
    vpa(my_sim.equations(:,3),4)
    
    % hardcoded from mocat3_eq_check
    disp('1.409e-5*D1^2 + 1.032e-6*D1*N1 + 2.818e-6*D1*S1 + 4.794e-8*N1^2 + 2.063e-7*N1*S1 - 221.3*N1 + 1.409e-7*S1^2 + 9.197*N2')
    disp('1.351e-5*D2^2 + 9.895e-7*D2*N2 + 2.703e-6*D2*S2 + 4.598e-8*N2^2 + 1.979e-7*N2*S2 - 9.197*N2 + 1.351e-7*S2^2 + 0.7843*N3')
    disp('1.297e-5*D3^2 + 9.499e-7*D3*N3 + 2.595e-6*D3*S3 + 4.414e-8*N3^2 + 1.9e-7*N3*S3 - 0.7843*N3 + 1.297e-7*S3^2 + 0.09009*N4')
    disp('1.246e-5*D4^2 + 9.127e-7*D4*N4 + 2.493e-6*D4*S4 + 4.241e-8*N4^2 + 1.825e-7*N4*S4 - 0.09009*N4 + 1.246e-7*S4^2 + 0.01518*N5')
    disp('1.198e-5*D5^2 + 8.775e-7*D5*N5 + 2.397e-6*D5*S5 + 4.078e-8*N5^2 + 1.755e-7*N5*S5 - 0.01518*N5 + 1.198e-7*S5^2')

    else

    % Su
    disp('Su')
    vpa(my_sim.equations(:,1),4)

    % hardcoded from mocat4_eq_check
    disp('lam_Su1 - 0.2*Su1 - 8.891e-8*N1*Su1 - 2.831e-11*Su1^2 - 2.831e-7*D1*Su1 - (2.831e-11*S1*Su1^2)/(S1 + Su1)')
    disp('lam_Su2 - 0.2*Su2 - 8.528e-8*N2*Su2 - 2.715e-11*Su2^2 - 2.716e-7*D2*Su2 - (2.715e-11*S2*Su2^2)/(S2 + Su2)')
    disp('lam_Su3 - 0.2*Su3 - 8.187e-8*N3*Su3 - 2.606e-11*Su3^2 - 2.607e-7*D3*Su3 - (2.606e-11*S3*Su3^2)/(S3 + Su3)')
    disp('lam_Su4 - 0.2*Su4 - 7.866e-8*N4*Su4 - 2.504e-11*Su4^2 - 2.505e-7*D4*Su4 - (2.504e-11*S4*Su4^2)/(S4 + Su4)')
    disp('lam_Su5 - 0.2*Su5 - 7.563e-8*N5*Su5 - 2.408e-11*Su5^2 - 2.408e-7*D5*Su5 - (2.408e-11*S5*Su5^2)/(S5 + Su5)')
    

    % S
    disp('S')
    vpa(my_sim.equations(:,2),4)
    
    % hardcoded from mocat4_eq_check
    disp('lam_S1 - 0.2*S1 - 8.891e-8*N1*S1 - 2.831e-7*D1*S1 - (2.831e-11*S1^2*Su1)/(S1 + Su1)')
    disp('lam_S2 - 0.2*S2 - 8.528e-8*N2*S2 - 2.716e-7*D2*S2 - (2.715e-11*S2^2*Su2)/(S2 + Su2)')
    disp('lam_S3 - 0.2*S3 - 8.187e-8*N3*S3 - 2.607e-7*D3*S3 - (2.606e-11*S3^2*Su3)/(S3 + Su3)')
    disp('lam_S4 - 0.2*S4 - 7.866e-8*N4*S4 - 2.505e-7*D4*S4 - (2.504e-11*S4^2*Su4)/(S4 + Su4)')
    disp('lam_S5 - 0.2*S5 - 7.563e-8*N5*S5 - 2.408e-7*D5*S5 - (2.408e-11*S5^2*Su5)/(S5 + Su5)')

   
    % D
    disp('D')
    vpa(my_sim.equations(:,3),4)
    
    % hardcoded from mocat4_eq_check
    disp('2.298*D2 - 55.29*D1 + 0.02*S1 + 0.02*Su1 + 8.889e-8*N1*S1 + 8.889e-8*N1*Su1 - 2.831e-8*D1^2 - 8.889e-9*D1*N1 + 2.831e-7*D1*S1 + 2.831e-7*D1*Su1')
    disp('0.1959*D3 - 2.298*D2 + 0.02*S2 + 0.02*Su2 + 8.527e-8*N2*S2 + 8.527e-8*N2*Su2 - 2.715e-8*D2^2 - 8.527e-9*D2*N2 + 2.715e-7*D2*S2 + 2.715e-7*D2*Su2')
    disp('0.02251*D4 - 0.1959*D3 + 0.02*S3 + 0.02*Su3 + 8.185e-8*N3*S3 + 8.185e-8*N3*Su3 - 2.606e-8*D3^2 - 8.185e-9*D3*N3 + 2.606e-7*D3*S3 + 2.606e-7*D3*Su3')
    disp('0.003794*D5 - 0.02251*D4 + 0.02*S4 + 0.02*Su4 + 7.864e-8*N4*S4 + 7.864e-8*N4*Su4 - 2.504e-8*D4^2 - 7.864e-9*D4*N4 + 2.504e-7*D4*S4 + 2.504e-7*D4*Su4')
    disp('0.02*S5 - 0.003794*D5 + 0.02*Su5 + 7.562e-8*N5*S5 + 7.562e-8*N5*Su5 - 2.408e-8*D5^2 - 7.562e-9*D5*N5 + 2.408e-7*D5*S5 + 2.408e-7*D5*Su5')
    
    
    % N
    disp('N')
    vpa(my_sim.equations(:,4),4)
    
    % hardcoded from mocat4_eq_check
    disp('9.197*N2 - 221.3*N1 + 2.063e-9*N1*S1 + 2.063e-9*N1*Su1 + 1.409e-8*S1*Su1 + 1.409e-5*D1^2 + 4.794e-8*N1^2 + 1.409e-8*Su1^2 + 1.032e-6*D1*N1 + 2.818e-8*D1*S1 + 2.818e-8*D1*Su1')
    disp('0.7843*N3 - 9.197*N2 + 1.979e-9*N2*S2 + 1.979e-9*N2*Su2 + 1.351e-8*S2*Su2 + 1.351e-5*D2^2 + 4.598e-8*N2^2 + 1.351e-8*Su2^2 + 9.895e-7*D2*N2 + 2.703e-8*D2*S2 + 2.703e-8*D2*Su2')
    disp('0.09009*N4 - 0.7843*N3 + 1.9e-9*N3*S3 + 1.9e-9*N3*Su3 + 1.297e-8*S3*Su3 + 1.297e-5*D3^2 + 4.414e-8*N3^2 + 1.297e-8*Su3^2 + 9.499e-7*D3*N3 + 2.595e-8*D3*S3 + 2.595e-8*D3*Su3')
    disp('0.01519*N5 - 0.09009*N4 + 1.825e-9*N4*S4 + 1.825e-9*N4*Su4 + 1.246e-8*S4*Su4 + 1.246e-5*D4^2 + 4.241e-8*N4^2 + 1.246e-8*Su4^2 + 9.126e-7*D4*N4 + 2.493e-8*D4*S4 + 2.493e-8*D4*Su4')
    disp('1.755e-9*N5*S5 - 0.01519*N5 + 1.755e-9*N5*Su5 + 1.198e-8*S5*Su5 + 1.198e-5*D5^2 + 4.078e-8*N5^2 + 1.198e-8*Su5^2 + 8.775e-7*D5*N5 + 2.397e-8*D5*S5 + 2.397e-8*D5*Su5')


    end
end    
