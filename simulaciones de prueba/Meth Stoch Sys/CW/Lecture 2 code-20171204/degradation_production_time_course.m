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
M=100;

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
    
%Calculate the mean over all 100 realisations
mean_A=mean(rec_vector,2);

%Calculate the analytically determined mean from the master equation
analytical_mean=k2nu*(1-exp(-k1*time_vec))/k1;

%Plot the first 10 trajectories
if M>=10
    plot(time_vec,rec_vector(:,1:10));
else
    plot(time_vec,rec_vector);
end

hold on
%Also plot on the mean over 100 realisations
[h1]=plot(time_vec,mean_A,'r','linewidth',5);

%Also plot the analytically determined mean
[h2]=plot(time_vec,analytical_mean,'k--','linewidth',5);

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of particles')

legend([h1,h2],['Mean of M=',num2str(M),' simulations'],'Analytical mean');

exportfig(gcf,...
            ['production_degradation_realisations_and_mean.eps'],...
            'Format','eps2',...
            'Width','20',...
            'Color','cmyk',...
            'Resolution',300,...
            'FontMode','fixed',...
            'FontSize',21);
%             'LineWidth',6);
        %Save as a .fig as well
        saveas(gcf,['production_degradation_realisations_and_mean.fig'],'fig');
        
        %Plot the figures from lecture 3 using the same code
%Make a new figure
figure

%Also plot the analytically determined mean
plot(time_vec,analytical_mean,'k--','linewidth',5);

%Set the axis size
axis([0 T_final 0 max(max(rec_vector(:,1),20))]);

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of particles')
exportfig(gcf,...
    ['production_degradation_mean.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',21);
        %Save as a .fig as well
        saveas(gcf,['production_degradation_mean.fig'],'fig');

%Make a new figure
figure

%plot the analytically determined mean
[h2]=plot(time_vec,analytical_mean,'k--','linewidth',5);

legend([h2],'Analytical mean');

hold on
%Plot the first trajectory
plot(time_vec,rec_vector(:,1));

%Set the axis size
axis([0 T_final 0 max(max(rec_vector(:,1),20))]);

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of particles')
exportfig(gcf,...
    ['production_degradation_mean_and_trajectory.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',21);
        %Save as a .fig as well
        saveas(gcf,['production_degradation_mean_and_trajectory.fig'],'fig');
