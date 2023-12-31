%--------------------------------------------------------------------------
% ETminUT: Difference ET-UT of ephemeris time and universal time
%
% Notes: The approximation spans the years from 1620 to 2050
%
% Input:
% x     the target year for interpolation
% Output:
% yp    interpolated value [s]
% 
% References:
% Meeus J.; Astronomical Algorithms; Willmann-Bell; Richmond, Virginia;
% 2nd edition (1998).
% https://eclipse.gsfc.nasa.gov/LEcat5/deltatpoly.html
%
% Last modified:   2018/01/27   Meysam Mahooti
%--------------------------------------------------------------------------
function yp = ETminUT(x)

if (2007 <= x && x <= 2050)
    t = x - 2000;
    yp = 62.92 + 0.32217 * t + 0.005589 * t^2;
    return
end

x1 = 1620;
xn = 2010;
y  = [...
	  121.0, 112.0, 103.0, 95.0, 88.0, 82.0, 77.0, 72.0, 68.0, 63.0,...
	  60.0, 56.0, 53.0, 51.0, 48.0, 46.0, 44.0, 42.0, 40.0, 38.0,...
	  35.0, 33.0, 31.0, 29.0, 26.0, 24.0, 22.0, 20.0, 18.0, 16.0,...
	  14.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 7.0, 7.0, 7.0,...
	  7.0, 7.0, 8.0, 8.0, 9.0, 9.0, 9.0, 9.0, 9.0, 10.0,...
	  10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 11.0, 11.0, 11.0,...
	  11.0, 11.0, 12.0, 12.0, 12.0, 12.0, 13.0, 13.0, 13.0, 14.0,...
	  14.0, 14.0, 14.0, 15.0, 15.0, 15.0, 15.0, 15.0, 16.0, 16.0,...
	  16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 15.0, 15.0, 14.0, 13.0,...
	  13.1, 12.5, 12.2, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 11.9,...
	  11.6, 11.0, 10.2, 9.2, 8.2, 7.1, 6.2, 5.6, 5.4, 5.3,...
	  5.4, 5.6, 5.9, 6.2, 6.5, 6.8, 7.1, 7.3, 7.5, 7.6,...
	  7.7, 7.3, 6.2, 5.2, 2.7, 1.4, -1.2, -2.8, -3.8, -4.8,...
	  -5.5, -5.3, -5.6, -5.7, -5.9, -6.0, -6.3, -6.5, -6.2, -4.7,...
	  -2.8, -0.1, 2.6, 5.3, 7.7, 10.4, 13.3, 16.0, 18.2, 20.2,...
	  21.1, 22.4, 23.5, 23.8, 24.3, 24.0, 23.9, 23.9, 23.7, 24.0,...
	  24.3, 25.3, 26.2, 27.3, 28.2, 29.1, 30.0, 30.7, 31.4, 32.2,...
	  33.1, 34.0, 35.0, 36.5, 38.3, 40.2, 42.2, 44.5, 46.5, 48.5,...
	  50.5, 52.2, 53.8, 54.9, 55.8, 56.9, 58.3, 60.0, 61.6, 63.0,...
	  63.8, 64.3, 64.6, 64.8, 65.5, 66.1];

if (length(y) > 3)
	interval = (xn - x1) / (length(y)-1);
	nearestX = fix((x-x1)/interval + .5);
	if (nearestX < 1)
		nearestX = 1;
    elseif (nearestX > length(y)-2)
		nearestX = length(y) - 2;
	end
	y = y(nearestX : nearestX+3);
	xn = x1 + (nearestX+1)*interval;
	x1 = x1 + (nearestX-1)*interval;
end

a = y(2) - y(1);
b = y(3) - y(2);
c = b - a;
abSum = a + b;
xSum = xn + x1;
xDiff = xn - x1;
n = (2*x - xSum) / xDiff;
yp = y(2) + n*.5*(abSum+n*c);

