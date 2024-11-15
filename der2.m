function [xx, yy] = der2(x, y, n)

if length(x) ~= length(y)
    error('x and y must be the same length')
end

if n < 3
    error('impossibile fra')
end

if nargin == 2
    n = 3;
end

l = length(x);
yy = zeros(1, l);

for i=1:l-n+1
    xv = x(i:i+n-1);
    yv = y(i:i+n-1);
    p = polyfit(xv, yv, 2);
    yy(i) = 2*p(1);
end

xx = x;