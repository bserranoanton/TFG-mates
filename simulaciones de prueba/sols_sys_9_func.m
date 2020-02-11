function [cSol,aSol,pSol,dSol] = sols_sys_9_func(lambda_taup,lambda_pp, r_tau, p0, lambda_pd, d0, mu_pc, c0, mu_da, a0)
    syms p(t)
    syms d(t)
    syms c(t)
    syms a(t)

    %Solve p:
    ode = diff(p,t) == lambda_taup*r_tau - lambda_pp*p;

    cond = p(0) == p0;
    pSol(t) = dsolve(ode,cond);

    %Solve d:
    ode = diff(d,t) == lambda_pd*pSol;

    cond = d(0) == d0;
    dSol(t) = dsolve(ode,cond);

    %Solve c:
    ode = diff(c,t) == mu_pc*pSol;

    cond = c(0) == c0;
    cSol(t) = dsolve(ode,cond);
    
    %Solve a:
    ode = diff(a,t) == mu_da*dSol;

    cond = a(0) == a0;
    aSol(t) = dsolve(ode,cond);
end