
%solucion del sistema (9)
syms c(t) a(t) p(t) d(t)

ode1 = diff(c) == -p;
ode2 = diff(a) == -d;
ode3 = diff(p) == -p;
ode4 = diff(d) == 0.5*p;
odes = [ode1; ode2;ode3;ode4];

cond1 = c(0) == 1;
cond2 = a(0) == 1;
cond3 = p(0) == 0.5;
cond4 = d(0) == 0;
conds = [cond1; cond2;cond3;cond4];

[cSol(t), aSol(t),pSol(t),dSol(t)] = dsolve(odes,conds);

disp(cSol);

fplot(cSol)
hold on
fplot(aSol)
hold on
fplot(pSol)
hold on
fplot(dSol)
grid on
legend('cSol','aSol','pSol','dSol','Location','best')