syms t_cell(t)
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
ode = diff(p,t) == alpha*p - beta*t_cell_value*p;

cond = p(0) == p0;
pSol(t) = dsolve(ode,cond);

fprintf('La solucion de p es: ');
disp(pSol);

%Resolvemos T:
ode = diff(t_cell,t,2) == -k*t_cell + lambda*p_value;
Dt_cell = diff(t_cell);
       
cond1 = t_cell(0) == t_cell0;
cond2 = Dt_cell(0) == t_cell_tilde0;

conds = [cond1 cond2];
t_cell_Sol(t) = dsolve(ode,conds);
t_cell_Sol(t) = simplify(t_cell_Sol);

fprintf('La solucion de T es: ');
disp(t_cell_Sol);