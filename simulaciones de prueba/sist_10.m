function dX=sist_10(t,X, lambda_pd)

dX(1)=-X(3);
dX(2)=-X(4);
dX(3)=-X(3);
dX(4)=lambda_pd*X(3);
dX=dX';