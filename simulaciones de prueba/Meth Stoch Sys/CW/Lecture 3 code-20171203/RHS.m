function dX=RHS(t,X,k1onu,k2onu,k3nu,k4nu)

dX(1)=-2*k1onu*X(1)^2-k2onu*X(1)*X(2)+k3nu;
dX(2)=-k2onu*X(1)*X(2)+k4nu;
dX=dX';