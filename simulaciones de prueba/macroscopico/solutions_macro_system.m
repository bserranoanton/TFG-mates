syms t_cell(t)
syms Dt_cell(t)
syms lambda
syms k
syms t_cell0
syms t_cell_tilde0
syms p_value

syms p(t)
syms alpha
syms beta
syms p0
syms t_cell_value

%Resolvemos p:
odeP = diff(p,t) == alpha*p - beta*t_cell*p;
condP = p(0) == p0;

%Resolvemos T:
odeT = diff(t_cell,t,2) == -k*t_cell + lambda*p;
Dt_cell = diff(t_cell);
       
condT1 = t_cell(0) == t_cell0;
condT2 = Dt_cell(0) == t_cell_tilde0;

eqns = [odeP, odeT];
conds = [condP, condT1, condT2];

sols(t) = dsolve(eqns,conds);

fprintf('La solucion de P es: ');
disp(sols.p);

fprintf('La solucion de T es: ');
disp(sols.t_cell);
