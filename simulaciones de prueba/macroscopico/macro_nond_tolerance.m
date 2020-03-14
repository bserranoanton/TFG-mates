syms t_cell(t) p(t) 
%Clohessy-Wiltshire Equations
% d2x = 2*n*dy + 3*(n^2)*x;
% d2y = -2*n*dx;
% d2z = (-n^2)*z;
%Constants
% a = 1.5;
% b = 0.1;
% 
% k = 0.4;
% lambda = 0.5;
% 
% a_star = a/sqrt(k);
% b_star = b*lambda/k*sqrt(k);

a_star = 1.1;
b_star = 0.01;

t0 = 0; %seconds
tf = 9.5; %seconds
dt_cell=diff(t_cell,t);

%Initial Conditions
c1 = 1; %p(0)
c2 = 0;  %T(0)
c3 = 0;  %T'(0)
y0 = [c1 c2 c3];
eq1 = diff(t_cell,t,2) == -t_cell + p;
eq2 = diff(p,t) == a_star*p - b_star*t_cell*p;

vars = [t_cell(t); p(t)]
[V,S] = odeToVectorField([eq1,eq2])


M = matlabFunction(V,'vars', {'t','Y'});
interval = [t0 tf];  %time interval    
% Impose a nonnegativity constraint 
% option1 = odeset('NonNegative',1);
 option2 = odeset('NonNegative',2); %T >= 0
% option3 = odeset('NonNegative',3);
% options = [option1, option2,option3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
ySol = ode45(M,interval,y0, option2);
tValues = linspace(interval(1),interval(2),1000);
yValues = deval(ySol,tValues,1); %number 1 denotes first position: pathogen

%plot 
figure 
xlabel('tiempo');  ylabel('Número de células');
[hA2] = plot(tValues,yValues,'r','LineWidth', 1);  %patógeno
% 
hold on
yValues = deval(ySol,tValues,2); 
[hA1] = plot(tValues,yValues,'b','LineWidth', 1);  %células T 
% ylim([0,1]);
%xlim([0,tf]);

%legend([hA2,hA1],'Patógeno','Células T');

