syms p(t)
syms lambda_taup_mem
syms lambda_pp_mem
syms r_tau
syms p0

syms c(t)
syms mu_pc_mem
syms c0

%Resolvemos p:
ode = diff(p,t) == lambda_taup_mem*r_tau - lambda_pp_mem*p;

cond = p(0) == p0;
pSol(t) = dsolve(ode,cond);

fprintf('La solucion de p es: ');
disp(pSol);

%Resolvemos c:
ode = diff(c,t) == -mu_pc_mem*pSol;

cond = c(0) == c0;
cSol(t) = dsolve(ode,cond);

fprintf('La solucion de c es: ');
disp(cSol);
