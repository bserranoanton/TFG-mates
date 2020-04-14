%[T, Y] = ode45(@YourFun, T, Y0, Opt);
function [position, isterminal, direction] = myEvent(t,ySol)
    position   = ySol(2);
    isterminal = 1;   % Stop the integration
    direction  = 0;
end