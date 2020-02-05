%This code is desgined to simulate stochatic production and degradation using an
%event-based algorithm
%By Kit Yates
%Created 25/08/15
%Last Modified 25/08/15

clear all
close all

%Define the rate constant of degradation
k1=0.1; %(sec^{-1})
k2nu=1; %(sec^{-1})

% Define the number of repeats we will do
M=10;

%Define A_e to the bt he point from which we exit the simulation
A_e=19;

%Define the vector of the sum of the first passage time over M repeats
time_sum=zeros(A_e+1,1);

%Define the initial number of particles
A_init=0;

%Define the outer for loop of all the different values of X_init we will run
%through
for A_init=0:A_e;
    %output X_init to see where we are
    A_init
    
    
    %Enter the inner for loop which runs through 1000 simulations for each
    %initial condition
    for jj=1:M
    
        %Define the initial time to be zero
        t=0;
        
        %Initialise X
        A=A_init;
        
        if A<A_e
            %Define loop_condition to be one initially so we enter the loop
            loop_condition=1;
        else
            %Provided we are not already outside the domain
            loop_condition=0;
        end
        
        %Run through a for loop for each of the time-steps
    %Enter the while loop
    while loop_condition
        
        %Calculate the propensity functions
        a0=A*k1+k2nu;
        
        %Determine the time for the next reaction
        tau=(1/a0)*log(1/rand);
        
        %Update the time
        t=t+tau;
        
        %Determine which reaction happened
        if rand*a0<k2nu
            %Implement the production reaction
            A=A+1;
        else
            %Implement the degradation reaction 
            A=max(A-1,0);
        end
    
         %If we have left then set the loop_condition to 0
            if A==A_e
                loop_condition=0;
            end
        
    end %while
        
        %Increment the time_Sum
        time_sum(A_init+1)=time_sum(A_init+1)+t;
        
    end %For jj
    
    
end %For A_init

%Calculate the average first passage times
MFPT=time_sum/M;

%Define how finely spaced our grid of x points should be
delta_A=0.01;

%Define a vector which will hold more finely spaced values of X
A_cont=0:delta_A:A_e;

%Find the length of this vector
length_A_cont=length(A_cont);

%Find the unnormalised SPD
SPD_CFPE=exp(-2*A_cont+(4*k2nu/k1-1)*log(k1*A_cont+k2nu));

%Normalise the SPD
SPD_CFPE=SPD_CFPE/(sum(SPD_CFPE*delta_A));

%Calculate the integral from 0 to z of the SPD
inner_int=cumsum(SPD_CFPE)*delta_A;

%Define a vector to hold the analytical MFPT
analytical_MFPT=zeros(length_A_cont,1);

for jj=1:length_A_cont
    %Calculate the analytical MFPT
    analytical_MFPT(jj)=delta_A*sum(inner_int(jj:end)./(SPD_CFPE(jj:end).*((k1*A_cont(jj:end)+k2nu)/2)));
end

%Define a plotting variable for A
A_plot=0:A_e;

%Make a new figure
figure

%Plot the SSA first passage time
bar(A_plot,MFPT,1);

hold on

%Plot the analytically derived MFTP on top
plot(A_cont,analytical_MFPT,'r','linewidth',5);

%Find the current axis sizes
axis_sizes=axis;

%Make new axis sizes
axis([ 0 A_e axis_sizes(3) axis_sizes(4)])

%Give x and y labels
xlabel('y')
ylabel('\tau(y) (sec)')

%Make a legend
legend('stochastic simulation','CFPE integral formula','location','SouthWest');

%
exportfig(gcf,...
    ['CFPE_exit_time_comparison.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',21);
%             'LineWidth',6);
%Save as a .fig as well
saveas(gcf,['CFPE_exit_time_comparison.fig'],'fig');
