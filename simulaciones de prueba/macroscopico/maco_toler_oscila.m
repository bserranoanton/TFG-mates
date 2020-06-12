syms t_cell(t) p(t) 
%Clohessy-Wiltshire Equations
% d2x = 2*n*dy + 3*(n^2)*x;
% d2y = -2*n*dx;
% d2z = (-n^2)*z;
%Constants
a = 0.02;
b = 0.1;

k = 0.2;
lambda = 0.1;

t0 = 0; %seconds
tf = 60; %seconds
dt_cell=diff(t_cell,t);

%Initial Conditions
c1 = 3; %p(0)
c2 = 0;  %T(0)
c3 = 0;  %T'(0)
y0 = [c1 c2 c3];
eq1 = diff(t_cell,t,2) == -k*t_cell + lambda*p;
eq2 = diff(p,t) == a*p - b*t_cell*p;

vars = [t_cell(t); p(t)]
[V,S] = odeToVectorField([eq1,eq2])


M = matlabFunction(V,'vars', {'t','Y'});
interval = [t0 tf];  %time interval    
% Impose a nonnegativity constraint 
option2 = odeset('NonNegative',2); %T >= 0



ySol = ode45(M,interval,y0, option2);
tValues = linspace(interval(1),interval(2),1000);
yValues = deval(ySol,tValues,1); %number 1 denotes first position likewise you can mention 2 to 6 for the next 5 positions

%plot 
figure 
xlabel('tiempo');  ylabel('Número de células');
[hA2] = plot(tValues,yValues/max(yValues),'r','LineWidth', 1);  %patógeno

hold on
yValues = deval(ySol,tValues,2); 
[hA1] = plot(tValues,yValues/max(yValues),'b','LineWidth', 1);  %células T 
ylim([0,1]);
set(gca,'YTickLabel',[]); %Para que no salgan los números del eje
set(gca,'XTickLabel',[]);


legend([hA2,hA1],'Patógeno','Células T');

