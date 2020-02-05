function dX=RHS_normal(t,X,k1onu,k2,A_init)

%X(1) is the first moment and X(2) the second

dX(2) = 4*k1onu*( 2.*(X(1))^3 - 3.*X(2).*X(1) + 2.*X(2) - X(1)) + 4*k2*(A_init.*(X(1) + 1) - X(2) - X(1)); 

dX=dX';