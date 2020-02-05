%This code is desgined to simulate the stationary probability distribution of the one-dimensional SDE
%8.21
%By Kit Yates
%Created 07/11/15
%Last Modified 07/11/15

clear all
close all

% Define the number of repeats we will do
M=10000;

%Define the timestep for the simulations
delta_t=10^(-3);

%Define the rate constants
k1=10^(-3); %(min^{-1})
k2=0.75; %(min^{-1})
k3=165; %(min^{-1})
k4=10^4; %(min^{-1})
k5=200;

%Define X_u to be the point to get to before we exit
X_u=250;

%Define the vector of the sum of the first passage time over M repeats
time_sum=zeros(X_u+1,1);

%Define a plotting variable for X
X_plot=0:X_u;


%Define the outer for loop of all the different values of X_init we will run
%through
for X_init=0:X_u;
    %output X_init to see where we are
    X_init
    
    
    %Enter the inner for loop which runs through 1000 simulations for each
    %initial condition
    for jj=1:M
    
        %Define the initial time to be zero
        t=0;
        
        %Initialise X
        X=X_init;
        
        if X<X_u
            %Define loop_condition to be one initially so we enter the loop
            loop_condition=1;
        else
            %Provided we are not already outside the domain
            loop_condition=0;
        end
        
        %Run through a while loop until the trajectory croses
        while loop_condition
            
            %Update time
            t=t+delta_t;
            
            %Rememebe the previous value of X
            X_old=X;
            
            %Run through the for loop for the time setps
            %Increment the vector of state variables
            X=X+delta_t*(-k1*X^3+k2*X^2-k3*X+k4)+sqrt(delta_t)*k5*randn;
            
            %If we have left then set the loop_condition to 0
            if X>X_u
                
                loop_condition=0;
            else
                %Check the extra condition in case the path has crossed
                if exp(-(X-X_u)*(X_old-X_u)/(((k5^2)/2)*delta_t))>rand
                    loop_condition=0;
                end
            end
            
        end %While
        
        %Increment the time_Sum
        time_sum(X_init+1)=time_sum(X_init+1)+t;
        
    end %For jj
    
    
end %For X_init

%Calculate the average first passage times
MFPT=time_sum/M;

%Define the size of the X increment
delta_X=0.01;

%Define a variable for the extent of the SPD
X_int=0:delta_X:2*X_u;

%Calculate the leength of X_int
length_X_int=length(X_int);

%Calculate the mean exit time using the stationary probability
%distribution
p1=exp((-3*k1*X_int.^4+4*k2*X_int.^3-6*k3*X_int.^2+12*k4*X_int)/(6*k5^2));

%normalise the solution
p1=p1/(sum(p1)*delta_X);

%Calculate the integral from 0 to z of the SPD
inner_int=cumsum(p1)*delta_X;

%calculate the last index to sum to
last_index=round(X_u/delta_X)+1;

%Define a vector to hold the analytical MFPT
analytical_MFPT=zeros(last_index,1);


%Calculate the SPD
for jj=1:last_index
    %Calculate the analytical MFPT
    analytical_MFPT(jj)=delta_X*sum(inner_int(jj:last_index)./p1(jj:last_index))/((k5^2)/2);
end

%Make a new figure
figure

%Plot the SSA first passage time
bar(X_plot,MFPT,1);

hold on

%Plot the analytically derived MFTP on top
plot(X_int(1:last_index),analytical_MFPT,'r','linewidth',5);

%Find the current axis sizes
axis_sizes=axis;

%Make new axis sizes
axis([ 0 X_u axis_sizes(3) axis_sizes(4)])

%Give x and y labels
xlabel('y')
ylabel('\tau(y)')

%Make a legend
legend('stochastic simulation','integral formula','location','SouthWest');


%
exportfig(gcf,...
    ['exit_time_as_a_function_of_y.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',21);
%             'LineWidth',6);
%Save as a .fig as well
saveas(gcf,['exit_time_as_a_function_of_y.fig'],'fig');
