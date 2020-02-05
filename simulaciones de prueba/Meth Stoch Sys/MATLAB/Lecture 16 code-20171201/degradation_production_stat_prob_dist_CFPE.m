%This code is desgined to calculate the stationary probability distribution
%of the production and degradation example and compare summary stats to the
%Chemical Fokker-planck equation
%By Kit Yates
%Created 25/08/15
%Last Modified 25/08/15

clear all
close all

%Define the rate constant of degradation
k1=0.1; %(sec^{-1})
k2nu=1; %(sec^{-1})

% Define the number of repeats we will do
M=1000;

% Define the final time we will simulate to
T_final=100;

%Define the recording time interval
rec_step=0.01;

%Define a vector of the time_points
time_vec=[0:rec_step:T_final];

%Calculate the number of recording steps required
num_rec_steps=T_final/rec_step;

%Define the vector which will reocrd the number of particlesat each time
%point
rec_vector=zeros(num_rec_steps+1,M);

%Define the initial number of particles
A_init=0;

%Define the initial number of molecules forthe recording vector
rec_vector(1,:)=A_init*ones(1,M);

%Run through a for loop for each of the repeats
for i=1:M

    %Define the initial time to be zero
    t=0;
    
    %Initialise the times to help recording
    t_after=0;
    
%   initialise the number of particles for this repeat
    A=A_init;
    
    %Run through a for loop for each of the time-steps
    %Enter the while loop
    while t<T_final
        
        %Calculate the propensity functions
        a0=A*k1+k2nu;
        
        %Determine the time for the next reaction
        tau=(1/a0)*log(1/rand);
        
        %Update the time
        t=t+tau;
        
        %Write the current value of A to A_old
        A_old=A;
        
        %Determine which reaction happened
        if rand*a0<k2nu
            %Implement the production reaction
            A=A_old+1;
        else
            %Implement the degradation reaction 
            A=max(A_old-1,0);
        end
        
        %Calculate the times for recording
        t_before=t_after;
        t_after=t;
   
        %Calculate the indices of the time step before and the current time
        %step in terms of recording
        ind_before=ceil((t_before+eps)/rec_step);
        ind_after=min(floor(t_after/rec_step),num_rec_steps);
        
        %Find out how many time-steps to write to
        steps_to_write=ind_after-ind_before+1;
        
        if steps_to_write>0 && steps_to_write~=Inf
           rec_vector(ind_before+1:ind_after+1,i)=(A_old)*ones(steps_to_write,1); 
        end
    
    end

end
    
%Define the time after which we consider the process to be stationary
t_stat=40;

%Disregard the first t_stat seconds worth of data
rec_vector_cut=rec_vector(40/rec_step:end,:);

%Find the maximum data value
max_data=max(rec_vector_cut(:));

%extract histogram data from each of the repeats after 40 seconds when the
%process should have roughly reached steady state
[N,X]=hist(rec_vector_cut(:),[0:1:max_data]);

%Normalised the data to make it a PDF
N=N/sum(N);

% Plot the histogram of the stationary frequencies
[h1]=bar(X,N);

%correct the axes
axis([0 max_data+0.5 0 ceil(max(N)*100)/100])

%Use the analytical formula in order to calculate the SPD for the
%corresponding CFPE

%Define how finely spaced our grid of x points should be
delta_X=0.01;

%Defie a vector which will hold more finely spaced values of X
X_cont=0:delta_X:max_data;

%Find the unnormalised SPD
SPD_CFPE=exp(-2*X_cont+(4*k2nu/k1-1)*log(k1*X_cont+k2nu));

%Normalise the SPD
SPD_CFPE=SPD_CFPE/(sum(SPD_CFPE*delta_X));

hold on
%Also plot the analytically determined mean
[h2]=plot(X_cont,SPD_CFPE,'r','linewidth',5);

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('number of particles')
ylabel('stationary distribution')

legend([h1,h2],['Long-run SSA'],'CFPE');

exportfig(gcf,...
            ['CFPE_SPD_comparison.eps'],...
            'Format','eps2',...
            'Width','20',...
            'Color','cmyk',...
            'Resolution',300,...
            'FontMode','fixed',...
            'FontSize',21);
%             'LineWidth',6);
        %Save as a .fig as well
        saveas(gcf,['CFPE_SPD_comparison.fig'],'fig');