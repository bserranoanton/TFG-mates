%This code is desgined to simulate the stochatic dimerisation degradation
%and production reactions (4.1).
% This will reproduce figure 4.1 (a)
%By Kit Yates
%Created 27/09/15
%Last Modified 27/09/15

clear all
close all

%Define the rate constant of degradation
k1onu=0.005; %(sec^{-1})
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

%Define the vector which will reocrd the number of particles at each time
%point
rec_vector_A=zeros(num_rec_steps+1,M);


%Define the initial number of particles
A_init=0;

%Define the initial number of molecules forthe recording vector
rec_vector_A(1,:)=A_init*ones(1,M);

%Define a vector which will hold the propensity functions
a=zeros(2,1);

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
        
        %Define the propensity functions
        a(1)=A*(A-1)*k1onu;
        a(2)=k2nu;
        
        %Calculate the cumulative sum of a
        cumsuma=cumsum(a);
        
        %Calculate the sum of the propensity functions
        a0=cumsuma(end);
        
        %Determine the time for the next reaction
        tau=(1/a0)*log(1/rand);
        
        %Update the time
        t=t+tau;
        
        %Write the current value of A to A_old
        A_old=A;
        
        %multiply the random number by the sum of the propensities and
        %store so we can reuse
        ra0=rand*a0;
        
        %if we have not run out of molecules
        if tau<inf;
            %Determine which reaction happened
            if ra0<cumsuma(1)
                %Implement the degradation reaction
                A=A_old-2;
            else
                % Implement a production reaction
                A=A_old+1;
            end
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
            rec_vector_A(ind_before+1:ind_after+1,i)=(A_old)*ones(steps_to_write,1);
        end
        
    end
    
end

%Calculate the mean over all 100 realisations
mean_A=mean(rec_vector_A,2);

%Solve the ODEs for the deterministic system
[t_deterministic,A_deterministic]=ode15s(@(t,X)RHS_dimerisation_degradation(t,X,k1onu,k2nu),[0,T_final],[0]);

%Plot the first 5 trajectories for A
f1=figure;
if M>=5
%     Change to figure f1
    figure(f1)
    plot(time_vec,rec_vector_A(:,1:10));
    
else
    
%     Change to figure f1
    figure(f1)
    plot(time_vec,rec_vector_A);
    
end

%     Change to figure f1
figure(f1)
hold on
%Also plot on the mean over 100 realisations
[hA1]=plot(time_vec,mean_A,'r','linewidth',5);

%Also plot the analytically determined mean
[hA2]=plot(t_deterministic,A_deterministic,'k--','linewidth',5);

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of particles')

legend([hA1,hA2],['Mean of M=',num2str(M),' simulations'],'Solution of ODE');

%     Change to figure f1
figure(f1)
exportfig(gcf,...
    ['dimerisation_degradation_ODE_and_paths.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',21);
%             'LineWidth',6);
%Save as a .fig as well
saveas(gcf,['dimerisation_degradation_ODE_and_paths.fig'],'fig');