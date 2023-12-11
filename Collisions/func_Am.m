function out = func_Am(d,ObjClass)
% NASA's new breakup model of evolve 4.0
% Author: N.L.Johnson, P.H.Krisko, J.-C. Liou, P.D.Anz-Meador
% Eq(6) pg 1381
% Distribution function for spacecraft fragments with Lc(=d)>11cm 

numObj = numel(d);      
logds = log10(d);        % d in meters
amsms = nan(numObj,5);  % store alpha,mu1,sig1,mu2,sig2

% if strcmpi(ObjClass,'Rocket Body') ||...
%    strcmpi(ObjClass,'Rocket Mission Related Object') ||...
%    strcmpi(ObjClass,'Rocket Fragmentation Debris') ||...
%    strcmpi(ObjClass,'Rocket Debris')
if ObjClass>4.5 && ObjClass<8.5 %Rocket-body related
    for ind = 1:numObj
        logd = logds(ind);
        % alpha
        if logd <= -1.4
            alpha = 1;
        elseif (-1.4 < logd) && (logd < 0)
            alpha = 1 - 0.3571 * (logd + 1.4);
        else % >= 0
            alpha = 0.5;
        end
    
        % mu1
        if logd <= -0.5
            mu1 = -0.45;
        elseif (-0.5 < logd) && (logd < 0)
            mu1 = -0.45 - 0.9 * (logd + 0.5);
        else % >= 0
            mu1 = -0.9;
        end
    
        % sigma1
        sigma1 = 0.55;
    
        % mu2
        mu2 = -0.9;
    
        % sigma2
        if logd <= -1.0
            sigma2 = 0.28;
        elseif (-1 < logd) && (logd < 0.1)
            sigma2 = 0.28 - 0.1636 * (logd + 1);
        else % >= 0.1
            sigma2 = 0.1;
        end
        amsms(ind,:) = [alpha, mu1, sigma1, mu2, sigma2];
    end
else        % not rocket body
    for ind = 1:numObj
        logd = logds(ind);
        % alpha
        if logd <= -1.95
            alpha = 0;
        elseif (-1.95 < logd) && (logd < 0.55)
            alpha = 0.3 + 0.4 * (logd + 1.2);
        else % >= 0.55
            alpha = 1;
        end
    
        % mu1
        if logd <= -1.1
            mu1 = -0.6;
        elseif (-1.1 < logd) && (logd < 0)
            mu1 = -0.6 - 0.318 * (logd + 1.1);
        else % >= 0
            mu1 = -0.95;
        end
    
        % sigma1
        if logd <= -1.3
            sigma1 = 0.1;
        elseif (-1.3 < logd) && (logd < -0.3)
            sigma1 = 0.1 + 0.2 * (logd + 1.3);
        else % >= -0.3
            sigma1 = 0.3;
        end
    
        % mu2
        if logd <= -0.7
            mu2 = -1.2;
        elseif (-0.7 < logd) && (logd < -0.1)
            mu2 = -1.2 - 1.333 * (logd + 0.7);
        else % >= -0.1
            mu2 = -2.0;
        end
    
        % sigma2
        if logd <= -0.5
            sigma2 = 0.5;
        elseif (-0.5 < logd) && (logd < -0.3)
            sigma2 = 0.5 - (logd + 0.5);
        else % >= -0.3
            sigma2 = 0.3;
        end
        amsms(ind,:) = [alpha, mu1, sigma1, mu2, sigma2];
    end
end

% N1 = mu1+sigma1*randn;
N1 = amsms(:,2) + amsms(:,3) .* randn(numObj,1);
% N2 = mu2+sigma2*randn;
N2 = amsms(:,4) + amsms(:,5) .* randn(numObj,1);

% out = (10.^(alpha * N1 + (1-alpha) * N2)); % eq(3.40)
out = (10.^(amsms(:,1) .* N1 + (1-amsms(:,1)) .* N2));

% Unit, seems to be m^2/kg per Fig 6.

end