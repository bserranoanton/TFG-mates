function dX=RHS_mean_field(t,X,k1onu,k2,A_init)

dX(1)= -2*k1onu.*X(1).*(X(1)-1)+2*k2.*(A_init - X(1));

dX=dX';