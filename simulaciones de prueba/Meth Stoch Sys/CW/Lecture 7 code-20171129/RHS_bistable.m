function dX=RHS_bistable(t,X,k1onu2,k2onu,k3,k4nu)

dX(1)=-k1onu2*X(1)^3+k2onu*X(1)^2-k3*X(1)+k4nu;
dX=dX';