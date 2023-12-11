function p = densityexp(h)
%
% Uses table from pp. 537 (Table 8-4) in Vallado for simple exponential
% atmospheric density model.
%
% Exponential model, U. S. Standard Atmosphere (1976) and CIRA-72
%
% h is the height above ellipsoid in km.
% density p is returned in kg/km^3
%
%if (h >= 150)
%    warning('h is >= 160 km, density calc is broken!');
%end
% these numbers for density units of kg/m^3, conversion to kg/km^3 is at
% the end of the function.

% original by mike shoemaker
% modified by t. kelecy to use "switch" logic rather than "if-then" 12/20/12
% modified 12/26/2012 by t. kelecy to "vectorize" computations
%
ncol = size(h,2);
p    = zeros(1,ncol);
h    = max(h,0);  %  check for h >= 0
for ii=1:ncol
    switch 1
        case (h(ii) >= 0 && h(ii) < 25)
            h0 = 0; p0 = 1.225; H = 7.249;
        case (h(ii) >= 25 && h(ii) < 30)
            h0 = 25; p0 = 3.899e-2; H = 6.349;
        case (h(ii) >= 30 && h(ii) < 40)
            h0 = 30; p0 = 1.774e-2; H = 6.682;
        case (h(ii) >= 40 && h(ii) < 50)
            h0 = 40; p0 = 3.972e-3; H = 7.554;
        case (h(ii) >= 50 && h(ii) < 60)
            h0 = 50; p0 = 1.057e-3; H = 8.382;
        case (h(ii) >= 60 && h(ii) < 70)
            h0 = 60; p0 = 3.206e-4; H = 7.714;
        case (h(ii) >= 70 && h(ii) < 80)
            h0 = 70; p0 = 8.770e-5; H = 6.549;
        case (h(ii) >= 80 && h(ii) < 90)
            h0 = 80; p0 = 1.905e-5; H = 5.799;
        case (h(ii) >= 90 && h(ii) < 100)
            h0 = 90; p0 = 3.396e-6; H = 5.382;
        case (h(ii) >= 100 && h(ii) < 110)
            h0 = 100; p0 = 5.297e-7; H = 5.877;
        case (h(ii) >= 110 && h(ii) < 120)
            h0 = 110; p0 = 9.661e-8; H = 7.263;
        case (h(ii) >= 120 && h(ii) < 130)
            h0 = 120; p0 = 2.438e-8; H = 9.473;
        case (h(ii) >= 130 && h(ii) < 140)
            h0 = 130; p0 = 8.484e-9; H = 12.636;
        case (h(ii) >= 140 && h(ii) < 150)
            h0 = 140; p0 = 3.845e-9; H = 16.149;
        case (h(ii) >= 150 && h(ii) < 180)
            h0 = 150; p0 = 2.070e-9; H = 22.523;
        case (h(ii) >= 180 && h(ii) < 200)
            h0 = 180; p0 = 5.464e-10; H = 29.740;
        case (h(ii) >= 200 && h(ii) < 250)
            h0 = 200; p0 = 2.789e-10; H = 37.105;
        case (h(ii) >= 250 && h(ii) < 300)
            h0 = 250; p0 = 7.248e-11; H = 45.546;
        case (h(ii) >= 300 && h(ii) < 350)
            h0 = 300; p0 = 2.418e-11; H = 53.628;
        case (h(ii) >= 350 && h(ii) < 400)
            h0 = 350; p0 = 9.518e-12; H = 53.298;
        case (h(ii) >= 400 && h(ii) < 450)
            h0 = 400; p0 = 3.725e-12; H = 58.515;
        case (h(ii) >= 450 && h(ii) < 500)
            h0 = 450; p0 = 1.585e-12; H = 60.828;
        case (h(ii) >= 500 && h(ii) < 600)
            h0 = 500; p0 = 6.967e-13; H = 63.822;
        case (h(ii) >= 600 && h(ii) < 700)
            h0 = 600; p0 = 1.454e-13; H = 71.835;
        case (h(ii) >= 700 && h(ii) < 800)
            h0 = 700; p0 = 3.614e-14; H = 88.667;
        case (h(ii) >= 800 && h(ii) < 900)
            h0 = 800; p0 = 1.170e-14; H = 124.64;
        case (h(ii) >= 900 && h(ii) < 1000)
            h0 = 900; p0 = 5.245e-15; H = 181.05;
        case (h(ii) >= 1000 )
            h0 = 1000; p0 = 3.019e-15; H = 268.00;
        otherwise
            error('h = %f',h(ii));
    end  %  switch 1
        
    p(ii) = p0 * exp((h0-h(ii))/H);
    
end  %  for ii=1:ncol
    
% convert density to kg/km^3
p = p ;%* 1000^3;
    
end  %  end of function