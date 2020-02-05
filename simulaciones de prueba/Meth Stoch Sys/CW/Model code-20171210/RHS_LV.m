function dX=RHS_LV(t,X,k1,k2onu,k3)

dX(1)=X(1)*(k1-k2onu*X(2));
dX(2)=X(2)*(k2onu*X(1)-k3);
dX=dX';