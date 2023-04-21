function [x,y] = NACA_Airfoils(m,p,t,c,N)
%%
% Input: m (maximum cambered line), p (location of maximum cambered line),
%        t (thickness of the airfoil), c (chord length), N (number of
%        points)
% Output: x-coordinate of chord, y-coordinate of chord
%%

point = linspace(0,c,N);
m = m/100;
p = p/10;
t = t/100;
y_t = @(x) t*c*(0.2969*sqrt(x/c) -0.126*(x/c) -0.3516*(x/c)^2 +0.2843*(x/c)^3 -0.1036*(x/c)^4)/0.2;
y_c1 = @(x) m*x*(2*p - x/c)/p^2;
y_c2 = @(x) m*(c - x)*(1 + (x/c) - 2*p)/(1 - p)^2;
dydx1 = @(x) 2*m*(1 - x/(c*p))/p;
dydx2 = @(x) 2*m*(p - x/c)/(1 - p)^2;

for i = 1:N
    if 0 < point(i) && point(i) <= p*c
        y_c = y_c1(point(i));
        theta = atan(dydx1(point(i)));
    else
        y_c = y_c2(point(i));
        theta = atan(dydx2(point(i)));
    end
    
    x_U(i) = point(i) - y_t(point(i))*sin(theta);
    x_L(i) = point(i) + y_t(point(i))*sin(theta);
    y_U(i) = y_c + y_t(point(i))*cos(theta);
    y_L(i) = y_c - y_t(point(i))*cos(theta);

end
x = [fliplr(x_L),x_U(2:end)];
y = [fliplr(y_L),y_U(2:end)];
x(isnan(x)) = 0;
y(isnan(y)) = 0;

end