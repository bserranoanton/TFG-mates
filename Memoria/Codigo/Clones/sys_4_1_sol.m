function [c,a,p,d] = sys_4_1_sol(t, lambda_taup,lambda_pp, r_tau, p0, lambda_pd, d0, mu_pc, c0, mu_da, a0)
p=(lambda_taup*r_tau + exp(-lambda_pp*t)*(lambda_pp*p0 - lambda_taup*r_tau))/lambda_pp;
d=d0 + (lambda_pd*p0 - (lambda_pd*lambda_taup*r_tau)/lambda_pp)/lambda_pp - ...
    (exp(-lambda_pp*t)*(lambda_pd*p0 - (lambda_pd*lambda_taup*r_tau)/lambda_pp) - ...
    lambda_pd*lambda_taup*r_tau*t)/lambda_pp;
c=c0 - (mu_pc*p0 - (lambda_taup*mu_pc*r_tau)/lambda_pp)/lambda_pp + ...
    (exp(-lambda_pp*t)*(mu_pc*p0 - (lambda_taup*mu_pc*r_tau)/lambda_pp) - ...
    lambda_taup*mu_pc*r_tau*t)/lambda_pp;
a=a0 + (lambda_pd*mu_da*p0 - (lambda_pd*lambda_taup*mu_da*r_tau)/lambda_pp)/lambda_pp^2 ...
    - (t*(d0*mu_da*lambda_pp^2 + lambda_pd*mu_da*p0*lambda_pp - ...
    lambda_pd*lambda_taup*mu_da*r_tau) + exp(-lambda_pp*t)*(lambda_pd*mu_da*p0 ...
    - (lambda_pd*lambda_taup*mu_da*r_tau)/lambda_pp) + ...
    (lambda_pd*lambda_pp*lambda_taup*mu_da*r_tau*t^2)/2)/lambda_pp^2;

end