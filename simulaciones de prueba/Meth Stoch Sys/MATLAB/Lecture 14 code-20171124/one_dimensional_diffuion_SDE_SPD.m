%This code is desgined to simulate the one-dimensional diffusion SDE
%8.8's and find the probability distribution at t=1
%By Kit Yates
%Created 07/11/15
%Last Modified 07/11/15

clear all
close all

% Define the number of repeats we will do
M=100000;

%Define how many repeats of this we will do. We do blocks of repeats because
%of memory constraints
num_reps=10;

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
rec_vector_X=zeros(1,M*num_reps);

%Define the initial condition for X
X_init=0;

for i=1:num_reps
    %Increment the vector of state variables
    rec_vector_X((i-1)*M+1:M*i)=D*sqrt(delta_t)*sum(randn(num_rec_steps,M),1);
end

%Find the maximum data value
max_data=ceil(max(rec_vector_X));
min_data=-max_data;

%Define how many boxes to split the histogram into
hist_box_num=40;

%extract histogram data from each of the repeats after t_stat seconds when the
%process should have roughly reached steady state
[N,X]=hist(rec_vector_X,linspace(min_data,max_data,hist_box_num));

%Normalised the data to make it a PDF
N=N/(sum(N)*2*max_data/hist_box_num);


%Calculate the PDF at t=1 according to the FPE
p1=1/sqrt(2*pi*T_final)*exp(-(X.^2)/(2*T_final));

% Plot the histogram of the stationary frequencies
[h1]=bar(X,N);

hold on

%Plot the deterministic trajecotry for the mean of X
[h2]=plot(X,p1,'r','linewidth',5);

%correct the axes
axis([min_data max_data 0 ceil(max(N)*100)/100])

%make a legend for the figure
legend([h1,h2],['SSA'],'FPE');

%Set the x and y labels
xlabel('x')
ylabel('p(x,1)')

    %
    exportfig(gcf,...
        ['PDF_for_diffusion_SDE.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['PDF_for_diffusion_SDE.fig'],'fig');
    