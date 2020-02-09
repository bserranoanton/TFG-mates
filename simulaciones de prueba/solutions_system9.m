syms p(t)
syms lambda_taup
syms lambda_pp
syms r_tau
syms p0

syms d(t)
syms lambda_pd
syms d0

syms c(t)
syms mu_pc
syms c0

syms a(t)
syms mu_da
syms a0

%lambda_taup = 6*10^(-5);
%lambda_pp = 0.1;
%r_tau = 0.2;

%Resolvemos p:
ode = diff(p,t) == lambda_taup*r_tau - lambda_pp*p;

cond = p(0) == p0;
pSol(t) = dsolve(ode,cond);

fprintf('La solucion de p es: ');
disp(pSol);

%Resolvemos d:
ode = diff(d,t) == lambda_pd*pSol;

cond = d(0) == d0;
dSol(t) = dsolve(ode,cond);

fprintf('La solucion de d es: ');
disp(dSol);

%Resolvemos c:
ode = diff(c,t) == mu_pc*pSol;

cond = c(0) == c0;
cSol(t) = dsolve(ode,cond);

fprintf('La solucion de c es: ');
disp(cSol);

%Resolvemos a:
ode = diff(a,t) == mu_da*dSol;

cond = a(0) == a0;
aSol(t) = dsolve(ode,cond);

fprintf('La solucion de a es: ');
disp(aSol);
