%--------------------------------------------------------------------------
% sign: returns absolute value of a with sign of b
%--------------------------------------------------------------------------
function result = sign_(a, b)

if (b>=0)
    result = abs(a);
else
    result = -abs(a);
end