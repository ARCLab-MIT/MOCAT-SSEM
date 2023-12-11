%--------------------------------------------------------------------------
%
% Chebyshev approximation of 3-dimensional vectors
%
% Inputs:
%     N       Number of coefficients
%     Ta      Begin interval
%     Tb      End interval
%     Cx      Coefficients of Chebyshev polyomial (x-coordinate)
%     Cy      Coefficients of Chebyshev polyomial (y-coordinate)
%     Cz      Coefficients of Chebyshev polyomial (z-coordinate)
%
% Last modified:   2018/01/27   M. Mahooti
% 
%--------------------------------------------------------------------------
function ChebApp = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz)

% Check validity
if ( (t<Ta) || (Tb<t) )
    error('ERROR: Time out of range in Cheb3D::Value\n');
end

% Clenshaw algorithm
tau = (2*t-Ta-Tb)/(Tb-Ta);  

f1 = zeros(1,3);
f2 = zeros(1,3);

for i=N:-1:2
    old_f1 = f1;
    f1 = 2*tau*f1-f2+[Cx(i),Cy(i),Cz(i)];
    f2 = old_f1;
end

ChebApp = tau*f1-f2+[Cx(1),Cy(1),Cz(1)];

