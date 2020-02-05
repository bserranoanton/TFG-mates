%This code is designed to simulate the direction of locusts of two types.
%The locusts change direction in response to meeting two locusts going in
%the opposite direction with rate r and also change direction randomly with
%rate 1-r.
%By Kit Yates
%Created 13/10/15
%Last Modified 13/10/15

clear all
close all


%Define the total number of individuals
N=20;

% Define the final time we will simulate to
T_final=350;

%Define the factor by which we wish to increase the speed of the system
speed_factor=60; %Since T_final is given in minutes and rates in sec^-1

%Define the rate at which particles switch from one state to the other
%spontaneously
r1=0.0225*speed_factor;

r2=10*0.0453*speed_factor;

r3=0;%0.1664*speed_factor;

r4=0*speed_factor;

%Define the recording time interval
rec_step=0.01;

%Define a vector of the time_points
time_vec=[0:rec_step:T_final];

%Calculate the number of recording steps required
num_rec_steps=T_final/rec_step;

%We only need to record the number of X1 since X1+X2 is conserved
rec_vector_X1=zeros(num_rec_steps+1,1);

%Define the initial number of particles
X1=round(N/2);

%Define the number of left moving individuals initially
X2=N-X1;

%Define the initial number of molecules forthe recording vector
rec_vector_X1(1)=X1;

%Define a vector which will hold the propensity functions
a=zeros(6,1);

%Run through a for loop for each of the repeats
%Define the initial time to be zero
t=0;

%Initialise the times to help recording
t_after=0;    

%Run through a for loop for each of the time-steps
%Enter the while loop
while t<T_final
    
    
    alr=(r1*X1/N+r2*X1*X2/N^2+(r3*X1*X2^2)/N^3+r4*X1*X2^3/N^4);
    arl=(r1*X2/N+r2*X1*X2/N^2+(r3*X2*X1^2)/N^3+r4*X2*X1^3/N^4);
    
    %Find the sum of the propensity functions
    a0=arl+alr;
    
    %Determine which reaction will take place
    if (rand*a0)<alr
        %Then convert a right moving individual to a left moving individual
        X1=X1-1;
        X2=X2+1;
    else
        %Then convert a left moving indicidual to a right-moving individual
        X2=X2-1;
        X1=X1+1;
    end
    %Determine the time for the next reaction
    tau=(1/a0)*log(1/rand);
    
    %Update the time
    t=t+tau;
    
    %Write the current value of A to A_old
    X1_old=X1;
    
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
        rec_vector_X1(ind_before+1:ind_after+1)=(X1_old)*ones(steps_to_write,1);
    end
        
end

    %Plot the single trajectory of the system
    plot(time_vec,2*(rec_vector_X1)/N-1);      
    
    hold on
    
    % %Set the title
    % title(['t= ', num2str(t*rec_step_increment)])
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('z')
    
    if r3>0
        
        exportfig(gcf,...
            ['model_switching_behaviour.eps'],...
            'Format','eps2',...
            'Width','20',...
            'Color','cmyk',...
            'Resolution',300,...
            'FontMode','fixed',...
            'FontSize',21);
        %             'LineWidth',6);
        %Save as a .fig as well
        saveas(gcf,['model_switching_behaviour.fig'],'fig');
        
    else
        exportfig(gcf,...
            ['model_switching_behaviour_r3=0.eps'],...
            'Format','eps2',...
            'Width','20',...
            'Color','cmyk',...
            'Resolution',300,...
            'FontMode','fixed',...
            'FontSize',21);
        %             'LineWidth',6);
        %Save as a .fig as well
        saveas(gcf,['model_switching_behaviour_r3=0.fig'],'fig');
        
    end