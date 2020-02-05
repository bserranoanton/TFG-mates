%This code is desgined to simulate the stationary probability distribution of the one-dimensional SDE
%8.21
%By Kit Yates
%Created 07/11/15
%Last Modified 07/11/15

clear all
close all

% Define the number of repeats we will do
M=100;

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

%Enter the for loop for the number of time steps
    
        %Increment the vector of state variables
        
        %     End the repeat for loop
    
    
%Define the time after which we consider the process to be stationary
t_stat=100;

%Disregard the first t_stat seconds worth of data
rec_vector_cut=rec_vector_X(t_stat/delta_t:end,:);

%Find the maximum data value
max_data=ceil(max(rec_vector_cut(:))/100)*100;

%extract histogram data from each of the repeats after t_stat seconds when the
%process should have roughly reached steady state
[N,X]=hist(rec_vector_cut(:),[0:1:max_data]);

%Normalised the data to make it a PDF
N=N/sum(N);

%Calculate the stationary probability distribution according to the
%deterministic solution
p1=exp((-3*k1*X.^4+4*k2*X.^3-6*k3*X.^2+12*k4*X)/(6*k5^2));

%normalise the solution
p1=p1/sum(p1);

figure

% Plot the histogram of the stationary frequencies
[h1]=bar(X,N);

hold on

%plot the probability distribution
[h2]=plot(X,p1,'r','linewidth',5);

%correct the axes
axis([0 max_data+0.5 0 ceil(max(N)*100)/100])

%make a legend for the figure
legend([h1,h2],['SSA'],'FPE','location','North');

%Set the x and y labels
xlabel('x')
ylabel('p(x,1)')

%
exportfig(gcf,...
    ['SPD_for_schlogl_SDE.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',21);
%             'LineWidth',6);
%Save as a .fig as well
saveas(gcf,['SPD_for_schlogl_SDE.fig'],'fig');
