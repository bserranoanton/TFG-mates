%This code is desgined to simulate the stochatic dimerisation reactions
%(3.7) and (3.8)
%By Kit Yates
%Created 28/08/15
%Last Modified 28/08/15

clear all
close all

%Define the rate constant of degradation
k1=2.5; %(sec^{-1})
k2onu=2; %(sec^{-1})
k3onu2=0.003; %(sec^{-1}) %(sec^{-1})

% Define the number of repeats we will do
M=1;

% Define the final time we will simulate to
T_final=50;

%Define the recording time interval
rec_step=0.01;

%Define a vector of the time_points
time_vec=[0:rec_step:T_final];

%Calculate the number of recording steps required
num_rec_steps=T_final/rec_step;

%Define the vector which will reocrd the number of particlesat each time
%point
rec_vector_A=zeros(num_rec_steps+1,M);
rec_vector_B=zeros(num_rec_steps+1,M);


%Define the initial number of particles
A_init=25;
B_init=25;


%Define the initial number of molecules forthe recording vector
rec_vector_A(1,:)=A_init*ones(1,M);
rec_vector_B(1,:)=B_init*ones(1,M);

%Define a vector which will hold the propensity functions
a=zeros(6,1);

%Run through a for loop for each of the repeats
for i=1:M
    
    %Define the initial time to be zero
    t=0;
    
    %Initialise the times to help recording
    t_after=0;
    
    %   initialise the number of particles for this repeat
    A=A_init;
    B=B_init;
    
    %Run through a for loop for each of the time-steps
    %Enter the while loop
    while t<T_final
        
        %Define the propensity functions
        a(1)=B*k1;
        a(2)=A*B*k2onu;
        a(3)=A*(A-1)*B*k3onu2;
        a(4)=A*k1;
        a(5)= A*B*k2onu;
        a(6) = B*(B-1)*A*k3onu2;
        
        %Calculate the cumulative sum of a
        cumsuma=cumsum(a);
        
        %Calculate the propensity functions
        a0=cumsuma(end);
        
        %Determine the time for the next reaction
        tau=(1/a0)*log(1/rand);
        
        %Update the time
        t=t+tau;
        
        %Write the current value of A to A_old
        A_old=A;
        B_old=B;
        
        %multiply the random number by the sum of the propensities and
        %store so we can reuse
        ra0=rand*a0;
        
        %if we have not run out of molecules
        if tau<inf;
           %Determine which reaction happened
            if ra0<cumsuma(1)
                %Implement the production reaction
                A=A_old+1;
                B = B_old -1;
            elseif ra0<cumsuma(2)
                %Implement the degradation reaction
                A=A_old+1;
                B=B_old-1;
            elseif ra0<cumsuma(3)
                % Implement a production reaction
                A=A_old+1;
                B=B_old-1;
            elseif ra0<cumsuma(4)
                % Implement a production reaction
                A=A_old-1;
                B = B_old +1;
            elseif ra0<cumsuma(5)
                % Implement a production reaction
                A=A_old-1;
                B = B_old+1;
            else
                % Implement the other production reaction
                A = A_old -1;
                B=B_old+1;
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
           rec_vector_B(ind_before+1:ind_after+1,i)=(B_old)*ones(steps_to_write,1);
        end
        
    end
    
end

%Calculate the mean over all 100 realisations
mean_A=mean(rec_vector_A,2);
%mean_B=mean(rec_vector_B,2);

%Solve the ODEs for the deterministic system
%[t_deterministic,X_deterministic]=ode15s(@(t,X)RHS(t,X,k1onu,k2onu,k3nu,k4nu),[0,T_final],[0,0]);

%Write A and B from X
% A_deterministic=X_deterministic(:,1)';
% B_deterministic=X_deterministic(:,2)';

%Plot the first 10 trajectories for A and B
f1=figure;
%f2=figure;
if M>=10
%     Change to figure f1
    figure(f1)
    plot(time_vec,rec_vector_A(:,1:10));
    
%     Change to figure f2
    %figure(f2)
    %plot(time_vec,rec_vector_B(:,1:10));
else
    
%     Change to figure f1
    figure(f1)
    plot(time_vec,rec_vector_A);
    
%     Change to figure f2
    %figure(f2)
    %plot(time_vec,rec_vector_B);
end

%     Change to figure f1
figure(f1)
hold on
%Also plot on the mean over 100 realisations
%[hA1]=plot(time_vec,mean_A,'r','linewidth',5);

%Also plot the analytically determined mean
%[hA2]=plot(t_deterministic,A_deterministic,'k--','linewidth',5);

% %Set the title
 %title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of A particles')

%legend([hA1,hA2],['Mean of M=',num2str(M),' simulations'],'Solution of ODEs');


%     Change to figure f2
% figure(f2)
% hold on
% %Also plot on the mean over 100 realisations
% [hB1]=plot(time_vec,mean_B,'r','linewidth',5);
% 
% %Also plot the analytically determined mean
% %[hB2]=plot(t_deterministic,B_deterministic,'k--','linewidth',5);

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
% xlabel('time (sec)')
% ylabel('number of B particles')
% 
% legend([hB1,hB2],['Mean of M=',num2str(M),' simulations'],'Solution of ODEs');
% 

%     Change to figure f1
% figure(f1)
% exportfig(gcf,...
%     ['dimerisation_ODE_and_trajectories_A.eps'],...
%     'Format','eps2',...
%     'Width','20',...
%     'Color','cmyk',...
%     'Resolution',300,...
%     'FontMode','fixed',...
%     'FontSize',21);
% %             'LineWidth',6);
% %Save as a .fig as well
% saveas(gcf,['dimerisation_ODE_and_trajectories_A.fig'],'fig');

%     Change to figure f2
% figure(f2)
% exportfig(gcf,...
%     ['dimerisation_ODE_and_trajectories_B.eps'],...
%     'Format','eps2',...
%     'Width','20',...
%     'Color','cmyk',...
%     'Resolution',300,...
%     'FontMode','fixed',...
%     'FontSize',21);
% %             'LineWidth',6);
% %Save as a .fig as well
% saveas(gcf,['dimerisation_ODE_and_trajectories_B.fig'],'fig');
