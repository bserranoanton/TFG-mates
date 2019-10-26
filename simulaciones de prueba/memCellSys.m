function dcdt = memCellSys(t,c,lambda_pd)

dcdt = [-c(3); -c(4);-c(3); lambda_pd*c(3)];


% function dydt =  memCellSys(t,y)
% %VDP1  Evaluate the van der Pol ODEs for mu = 1
% %
% %   See also ODE113, ODE23, ODE45.
% 
% %   Jacek Kierzenka and Lawrence F. Shampine
% %   Copyright 1984-2014 The MathWorks, Inc.
% 
% dydt = [y(2); (1-y(1)^2)*y(2)-y(1)];