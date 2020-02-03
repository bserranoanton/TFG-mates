%Define the rate constant of degradation
k1onu2=4*10^(-5)*speed_factor; %(sec^{-1})
k2nu=50*speed_factor; %(sec^{-1})
k3=10*speed_factor; %(sec^{-1}
k4nu=25*speed_factor; %(sec^{-1})

% Define the final time we will simulate to
T_final=80;

%Define the initial number of particles
A_init=10;
B_init=10;

%Solve the ODEs for the deterministic system
%[t_deterministic,X_deterministic]=ode15s(@(t,X)RHS(t,X,k1onu2,k2nu,k3,k4nu),[0,T_final],[A_init,B_init]);
[t_deterministic,X_deterministic]=ode15s(@(t,X)RHS(t,X,k1onu2,k2nu,k3,k4nu),[0,T_final],[A_init,B_init]);
%Write A and B from X
A_deterministic=X_deterministic(:,1)';
B_deterministic=X_deterministic(:,2)';

figure
%plot the determinisitc trajectory fo A
[hA2]=semilogy(t_deterministic,A_deterministic,'k--','linewidth',5);
%Set the x and y labels
xlabel('time (min)')
ylabel('number of A particles')

%Find the original axis
original_axis=axis;

%Make sure the axes only go up to th final time
axis([0 80 original_axis(3) original_axis(4)])


%legend([hA1,hA2],'Stochastic trajectory','Deterministic trajectory');

%
exportfig(gcf,...
    ['resonance_stochastic_deterministic_A_log_scale.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',21);
%             'LineWidth',6);
%Save as a .fig as well
saveas(gcf,['resonance_stochastic_deterministic_A_log_scale.fig'],'fig');