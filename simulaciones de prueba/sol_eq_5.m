%function ySol = sol_eq_5(alpha, beta, N, y0)
  syms y(t)
  syms N
  syms y0
  syms beta
  syms alpha
  
  %Solve y:
  ode = diff(y,t) == (alpha - beta*N)*y;

  cond = y(0) == y0;
  ySol(t) = dsolve(ode,cond);
  
  disp(ySol);
%end