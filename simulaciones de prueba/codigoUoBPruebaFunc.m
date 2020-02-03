function dX=RHS(t,X,k1onu2,k2nu,k3,k4nu)

dX(1)=k1onu2*X(1)^2*X(2)+k2nu-k3*X(1);
dX(2)=-k1onu2*X(1)^2*X(2)+k4nu;
dX=dX';