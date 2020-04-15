%This code is desgined to simulate system 5.1.
%By Belen Serrano Anton
%Created 15/04/2020
%Last Modified 15/04/2020

function [max_P, max_T, t_max_P, t_max_T, t_min_P, t_min_T] = macro_func_neural_network(a, b, k, lambda)
syms t_cell(t) p(t) 

max_P = -1; max_T = -1;
t_max_P = -1; t_max_T = -1;
t_min_P = -1; t_min_T = -1;

t0 = 0; 
tf = 15; 
dt_cell=diff(t_cell,t);

%Initial Conditions
c1 = 3;     %P(0)
c2 = 1;     %T(0)
c3 = 0;     %T'(0)
y0 = [c1 c2 c3];
eq1 = diff(t_cell,t,2) == -k*t_cell + lambda*p;
eq2 = diff(p,t) == a*p - b*t_cell*p;

vars = [t_cell(t); p(t)];
[V,S] = odeToVectorField([eq1,eq2]);

M = matlabFunction(V,'vars', {'t','Y'});
interval = [t0 tf];  %Time interval    
%Impose a nonnegativity constraint 
option1 = odeset('Events', @myEvent);
%option2 = odeset('NonNegative',2); %T >= 0

[t,ySol,te,ye,ie] = ode45(M,interval,y0, option1);
tValues = linspace(interval(1),interval(2),1000);

%--------------------PATOGENO-----------------------------------------------
yValuesP = ySol(:,1);

%Hacemos que una vez que el patógeno es 0 no vuelva a reproducirse
flag_defeated = yValuesP < 0.01;
yValuesP(flag_defeated) = 0;

min_p_cell = min(yValuesP);
index_time_min_p_cell = yValuesP == min_p_cell;
time_min_p_cell = min(t(index_time_min_p_cell));

max_p_cell = max(yValuesP);
index_time_max_p_cell = yValuesP == max_p_cell;
time_max_p_cell = t(index_time_max_p_cell);

%--------------------CELULAS T----------------------------------------------

yValuesT = ySol(:,2);
min_t_cell = min(yValuesT);
index_time_min_t_cell = yValuesT == min_t_cell;
tValues(1) = 500; %Evita que el tiempo donde las células T son 0 sea 0
time_min_t_cell = min(t(index_time_min_t_cell));
tValues(1) = 0;
 
max_t_cell = max(yValuesT);
index_time_max_t_cell = yValuesT == max_t_cell;
time_max_t_cell = t(index_time_max_t_cell);


%If pathogen molecules are least than 0.01 we consider that pathogen has
%been totally defeated
if(yValuesP(index_time_min_t_cell) <= 0.01)%Intolerance
    max_P = max_p_cell;
    max_T = max_t_cell;
    t_max_P = time_max_p_cell;
    t_max_T = time_max_t_cell;
    t_min_P = time_min_p_cell;
    t_min_T = time_min_t_cell;
end



