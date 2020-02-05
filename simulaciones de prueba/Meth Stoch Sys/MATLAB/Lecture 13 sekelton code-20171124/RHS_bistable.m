function dX=RHS_bistable(t,X,k1,k2,k3,k4)

dX(1)=-k1*X(1)^3+k2*X(1)^2-k3*X(1)+k4;
dX=dX';