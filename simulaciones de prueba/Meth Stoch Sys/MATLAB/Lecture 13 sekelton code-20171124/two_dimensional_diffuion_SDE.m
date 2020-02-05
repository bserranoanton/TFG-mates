%This code is desgined to simulate the two-dimensional diffusion SDE
%8.18-8.19 in lecture 8
%By Kit Yates
%Created 05/011/15
%Last Modified 05/11/15

clear all
close all

% Define the number of repeats we will do
M=1000;

% Define the final time we will simulate to
T_final=10*60;

%Define the timestep for the simulations
delta_t=10^(-1);

%Define the diffusion coefficient
D=10^(-4);

%Define a vector of the time_points
time_vec=[0:delta_t:T_final];

%Define the  number of recording steps
num_rec_steps=T_final/delta_t;

%Define the vector which will reocrd the value of X at each time
%point
rec_vector_X=zeros(num_rec_steps+1,M);
rec_vector_Y=zeros(num_rec_steps+1,M);

%Define the initial condition for X
X_init=0;
Y_init=0;

%Write the initial condisiotn to this vector
rec_vector_X(1,:)=X_init*ones(1,M);
rec_vector_Y(1,:)=Y_init*ones(1,M);

    %Define the initial time to be zero
    t=0;
    
   
    
    %Increment the vector of state variables
    
    
    %Find the mean value of X and Y with time


    figure 
    
    %Plot the deterministic trajecotry for the mean of X
    [hXm]=plot(mean_X,mean_Y,'r','linewidth',5);
    
%     Finally plot 6 realisations of the stochastic trajectories
    [hXs]=plot(rec_vector_X(:,1:6),rec_vector_Y(:,1:6));
    
    hold on
       
    [hXf]=plot(rec_vector_X(end,1:6),rec_vector_Y(end,1:6),'ko','markersize',5,'linewidth',5)
    
    
    
    %Set the x and y labels
    xlabel('x (mm)')
    ylabel('y (mm)')

    %
    exportfig(gcf,...
        ['diffusive_trajectories_2_D.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['diffusive_trajectories_2_D.fig'],'fig');
    