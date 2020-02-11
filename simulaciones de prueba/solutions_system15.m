syms p(t)
syms lambda_taup
syms lambda_pp
syms r_tau
syms p0

syms c(t)
syms mu_pc
syms c0

%Resolvemos p:
ode = diff(p,t) == lambda_taup*r_tau - lambda_pp*p;

cond = p(0) == p0;
pSol(t) = dsolve(ode,cond);

fprintf('La solucion de p es: ');
disp(pSol);

%Resolvemos c:
ode = diff(c,t) == -mu_pc*pSol;

cond = c(0) == c0;
cSol(t) = dsolve(ode,cond);

fprintf('La solucion de c es: ');
disp(cSol);
