function [correct] = mocat4S_verify_equations()
% returns True if the hardcoded equations and workbook equations match
% returns False otherwise

N_shell = 5;

disp('Computing hardcoded equations ... ')
[hardcoded_eqs, hardcoded_results] = mocat4s_eq_check();
disp('Finished computing hardcoded equations.')


disp('Computing workbook equations ... ')
S = sym("S_",[N_shell, 1]);
Su = sym("Su_",[N_shell, 1]);
assume(S > 0);
assume(Su > 0);
workbook_sim = MOCAT_4S_Workbook_func();
disp('Finished computing workbook equations.')
%%
disp('Hardcoded equations:');
vpa(hardcoded_eqs, 3)

disp('Workbook equations:')
workbook_eqs = reshape(workbook_sim.equations, [20 1]);

vpa(workbook_eqs, 3)

disp('Equation Difference:')
vpa(hardcoded_eqs - workbook_eqs, 2)
diff = hardcoded_eqs - workbook_eqs;

% process equation so that coefficients can be extracted easily:
diff = subs(diff, S + Su, ones([N_shell, 1]));

tolerance = 1e-13;
correct = true;
for i=1:length(diff)-1
    [coefficients,t] = coeffs(diff(i,1));
    for c=1:length(coefficients)-1
        if coefficients(c) > tolerance 
            correct = false;
            return;
        end
    end
end

end