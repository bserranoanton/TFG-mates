function dX=RSH_dimerisation_degradation(t,X,k1onu,k2nu)

dX(1) = -2*k1onu*X(1)^2+k2nu;
dX = dX';
