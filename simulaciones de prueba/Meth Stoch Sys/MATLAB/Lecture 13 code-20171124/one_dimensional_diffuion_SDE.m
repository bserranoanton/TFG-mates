%This code is desgined to simulate the one-dimensional diffusion SDE
%8.8 in lecture 8
%By Kit Yates
%Created 05/011/15
%Last Modified 07/11/15

clear all
close all

% Define the number of repeats we will do
M=1000;

% Define the final time we will simulate to
T_final=1;

%Define the timestep for the simulations
delta_t=10^(-3);

%Define the diffusion coefficient
D=1;

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
    
   
    
    %Increment the vector of state variables
    rec_vector_X(2:end,:)=D*sqrt(delta_t)*cumsum(randn(num_rec_steps,M),1);
    
    %Find the mean value of X with time
    mean_X=mean(rec_vector_X,2);

    figure
    %Plot the mean of the stochastic trajectory for X
    [hXm]=plot(time_vec,mean_X,'r','linewidth',5);
    
    hold on
    
    %Plot the deterministic trajecotry for the mean of X
    [hXd]=plot(time_vec,zeros(1,num_rec_steps+1),'k--','linewidth',5);
    
%     Finally plot 6 realisations of the stochastic trajectories
    [hXs]=plot(time_vec,rec_vector_X(:,1:6));
    
    
    %make a legend for the figure
    legend([hXm,hXd],['Mean of M=',num2str(M),' simulations'],'Solution of ODE');
    
    %Set the x and y labels
    xlabel('time')
    ylabel('x')

    %
    exportfig(gcf,...
        ['diffusive_trajectories_1_D.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['diffusive_trajectories_1_D.fig'],'fig');
    