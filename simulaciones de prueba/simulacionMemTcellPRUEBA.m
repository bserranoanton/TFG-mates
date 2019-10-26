%Belén Serrano Antón

%Condiciones iniciales:
y0 = [1,1,1,1];

tspan = [-5 5];
[t,y] = ode45(@(t,y) memCellSys(t,y,0.5), tspan , y0);

plot(t,y(:,1),'-o',t,y(:,2),'-o',t,y(:,3),'-o',t,y(:,4),'-o');
title('Solution of van der Pol Equation (\mu = 1) with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2','y_3','y_4')


% [t,y] = ode45(@memCellSys,[0 20],[2; 0]);
% 
% plot(t,y(:,1),'-o',t,y(:,2),'-o')
% title('Solution of van der Pol Equation (\mu = 1) with ODE45');
% xlabel('Time t');
% ylabel('Solution y');
% legend('y_1','y_2')