syms x(t) y(t) z(t)
%Clohessy-Wiltshire Equations
% d2x = 2*n*dy + 3*(n^2)*x;
% d2y = -2*n*dx;
% d2z = (-n^2)*z;
%Constants
a = 6793.137; %km
mu = 398600.5; %km^3/s^2
n = sqrt(mu/a^3);
t0 = 0; %seconds
tf = 5400; %seconds
dx=diff(x,t);
dy=diff(y,t);
dz=diff(z,t);
%Initial Conditions
c1 = -0.066538073651029; %km
c2 =0.186268907590665; %km
c3 =0.000003725378152; %km
c4 = -0.000052436200437; %km/s
c5 =0.000154811363681; %km/s
c6 = 0.000210975508077; %km/s
y0 = [c1 c2 c3 c4 c5 c6];
eq1 = diff(x,2) == 2*n*dy + 3*(n^2)*x;
eq2 = diff(y,2) == -2*n*dx;
eq3 = diff(z,2) == (-n^2)*z;
vars = [x(t); y(t); z(t)]
V = odeToVectorField([eq1,eq2,eq3])
M = matlabFunction(V,'vars', {'t','Y'});
interval = [t0 tf];  %time interval    
ySol = ode45(M,interval,y0);
tValues = linspace(interval(1),interval(2),1000);
yValues = deval(ySol,tValues,1); %number 1 denotes first position likewise you can mention 2 to 6 for the next 5 positions
plot(tValues,yValues)   %if you want to plot position vs time just swap here