%Define the constant of ...
lambda_pd=0.5;
p0=2;
d0=1;
% Define the final time we will simulate to
T_final=10;

%Define the initial number of particles
c_init=1;
a_init=1;
p_init=p0;
d_init=d0;

%Solve the ODEs for the deterministic system
[t_deterministic,X_deterministic]=ode15s(@(t,X)sist_10(t,X,lambda_pd),[0,T_final],[c_init,a_init,p_init,d_init]);

%Write c,a,p,d from X
c_deterministic=X_deterministic(:,1)';
a_deterministic=X_deterministic(:,2)';
p_deterministic=X_deterministic(:,3)';
d_deterministic=X_deterministic(:,4)';

figure
%plot the determinisitc trajectory of c,a,p and d
%[hC2]=semilogy(t_deterministic,c_deterministic,'g','linewidth',5);
%[hA2]=semilogy(t_deterministic,a_deterministic,'b','linewidth',5);
%[hP2]=semilogy(t_deterministic,p_deterministic,'m','linewidth',5);
%[hD2]=semilogy(t_deterministic,d_deterministic,'c','linewidth',5);

[hC2] = plot(t_deterministic,c_deterministic,'g');
hold on
[hA2] = plot(t_deterministic,a_deterministic,'b');
hold on
[hP2] = plot(t_deterministic,p_deterministic,'m');
hold on
[hD2] = plot(t_deterministic,d_deterministic,'c');

%Set the x and y labels
xlabel('time (min)')
ylabel('number of particles')

%Find the original axis
original_axis=axis;

%Make sure the axes only go up to th final time
axis([0 10 original_axis(3) original_axis(4)])


legend([hC2,hA2,hP2,hD2],'c','a','p','d');

% %
% exportfig(gcf,...
%     ['resonance_stochastic_deterministic_A_log_scale.eps'],...
%     'Format','eps2',...
%     'Width','20',...
%     'Color','cmyk',...
%     'Resolution',300,...
%     'FontMode','fixed',...
%     'FontSize',21);
%             'LineWidth',6);
%Save as a .fig as well
%saveas(gcf,['resonance_stochastic_deterministic_A_log_scale.fig'],'fig');