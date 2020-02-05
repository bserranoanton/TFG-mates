%This code is desgined to simulate the one-dimensional diffusion SDE
%8.21 in lecture 8
%By Kit Yates
%Created 07/11/15
%Last Modified 07/11/15

clear all
close all

% Define the number of repeats we will do
M=200;

% Define the final time we will simulate to
T_final=800;

%Define the timestep for the simulations
delta_t=10^(-3);

%Define the rate constants
k1=10^(-3); %(min^{-1})
k2=0.75; %(min^{-1})
k3=165; %(min^{-1})
k4=10^4; %(min^{-1})
k5=200;

%Define a vector of the time_points
time_vec=[0:delta_t:T_final];

%Define the  number of recording steps
num_rec_steps=T_final/delta_t;

%Define the vector which will reocrd the value of X at each time
%point
rec_vector_X=zeros(num_rec_steps+1,M);

%Define the initial condition for X
X_init=0;

%Write the initial condisiotn to this vector
rec_vector_X(1,:)=X_init*ones(1,M);

    %Define the initial time to be zero
    t=0;
   
    for i=1:num_rec_steps
        %Increment the vector of state variables
        rec_vector_X(i+1,:)=rec_vector_X(i,:)+delta_t*(-k1*rec_vector_X(i,:).^3+k2*rec_vector_X(i,:).^2-k3*rec_vector_X(i,:)+k4)+sqrt(delta_t)*k5*randn(1,M);
    end
    
    %Find the mean value of X with time
    mean_X=mean(rec_vector_X,2);
    
    %Solve the ODEs for the deterministic system
    [t_deterministic,X_deterministic]=ode15s(@(t,X)RHS_bistable(t,X,k1,k2,k3,k4),[0,T_final],X_init);

    figure
    
   
    %Plot the deterministic trajecotry for the mean of X
    [hXd]=plot(t_deterministic,X_deterministic,'k--','linewidth',5);
   
     hold on
    
    %Plot the mean of the stochastic trajectory for X
    [hXm]=plot(time_vec,mean_X,'r','linewidth',5);
    
    
%     Finally plot 1 realisation of the stochastic trajectories
    [hXs]=plot(time_vec,rec_vector_X(:,1));
    
    %make a legend for the figure
    legend([hXm,hXd],['Mean of M=',num2str(M),' simulations'],'Solution of ODE','location','NorthWest');
    
    %Set the x and y labels
    xlabel('time')
    ylabel('x')

    %
    exportfig(gcf,...
        ['diffusion_and_drift_trajectories_1_D.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['diffusion_and_drift_trajectories_1_D.fig'],'fig');
    