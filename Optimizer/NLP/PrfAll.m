function [f,df] = PrfAll(x,data)

N_shell = data.N_shell;
S_all = x(1:N_shell);

% HIGH CAPACITY
S_all = S_all(1:14);
f = -sum(log(S_all));

% NOT HIGH CAPACITY
% f = - sum(log(S_all(6:8))); % shells at 500-600 km
% f = - sum(S_all(7)); % shell at 550 km

df = [];

end