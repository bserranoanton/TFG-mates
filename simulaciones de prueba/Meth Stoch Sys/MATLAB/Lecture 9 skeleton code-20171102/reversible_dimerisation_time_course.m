%This code is desgined to simulate the reversible dimerisation reactions
%(9.1)
%By Kit Yates
%Created 22/12/16
%Last Modified 22/12/16

clear all
close all

%Define the rate constant of degradation
k1onu=0.02/2; %(sec^{-1})
k2=0.2/2; %(sec^{-1})

% Define the number of repeats we will do
M=1000;

% Define the final time we will simulate to
T_final=25;

%Define the recording time interval
rec_step=0.01;

%Define a vector of the time_points
time_vec=[0:rec_step:T_final];

%Calculate the number of recording steps required
num_rec_steps=T_final/rec_step;

%Define the vector which will reocrd the number of particlesat each time
%point
rec_vector_A=zeros(num_rec_steps+1,M);

%Define the initial number of particles
A_init=50;

%Define the initial number of molecules for the recording vector
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
         a(2)=(A_init - A)*k2; %B/??
        
        %Calculate the sum of the propensity functions
        cumsuma = cumsum(a);
        
       
          %?????????Calculate the propensity functions
            a0=cumsuma(end);
        
        
        %Determine the time for the next reaction
        tau=(1/a0)*log(1/rand);

        
        %Update the time
        t=t+tau;
        
        %Write the current value of A to A_old
        A_old=A;
        
        
        
        %if we have not run out of molecules
        if tau<inf;
            %Determine which reaction happened
            if rand*a0<a(1)    %Implement the production reaction
                A=A_old-2;
            else
                A=A_old+2; %Implement the reverse dimerisation process
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

%Calculate the mean over all M realisations
mean_A=mean(rec_vector_A,2);

%Calculate the variance over all M realisations
var_A=var(rec_vector_A')';

%Solve the ODEs for the mean-field model system
[t_deterministic_mean_field,X_deterministic_mean_field]=ode15s(@(t,X)RHS_mean_field(t,X,k1onu,k2,A_init),[0,T_final],[A_init]);

%Write A from X
A_deterministic_mean_field=X_deterministic_mean_field(:,1)';

%Solve the ODEs for the pairwise model system
[t_deterministic_pairwise,X_deterministic_pairwise]=ode15s(@(t,X)RHS_pairwise(t,X,k1onu,k2,A_init),[0,T_final],[A_init,A_init^2]);

%Write A from X
A_deterministic_pairwise=X_deterministic_pairwise(:,1)';
A_deterministic_pairwise_second_moment=X_deterministic_pairwise(:,2)';

%Solve the ODEs for the normal closure system
[t_deterministic_normal,X_deterministic_normal]=ode15s(@(t,X)RHS_normal(t,X,k1onu,k2,A_init),[0,T_final],[A_init,A_init^2]);

%Write A from X
A_deterministic_normal=X_deterministic_normal(:,1)';
A_deterministic_normal_second_moment=X_deterministic_normal(:,2)';

%Solve the ODEs for the fourth-order closure system
[t_deterministic_fourth_cumulant,X_deterministic_fourth_cumulant]=ode15s(@(t,X)RHS_fourth_cumulant(t,X,k1onu,k2,A_init),[0,T_final],[A_init,A_init^2,A_init^3]);

%Write A from X
A_deterministic_fourth_cumulant=X_deterministic_fourth_cumulant(:,1)';
A_deterministic_fourth_cumulant_second_moment=X_deterministic_fourth_cumulant(:,2)';
A_deterministic_fourth_cumulant_third_moment=X_deterministic_fourth_cumulant(:,3)';

%Solve the ODEs for the kirkwood closure system
[t_deterministic_kirkwood,X_deterministic_kirkwood]=ode15s(@(t,X)RHS_kirkwood(t,X,k1onu,k2,A_init),[0,T_final],[A_init,A_init^2]);

%Write A from X
A_deterministic_kirkwood=X_deterministic_kirkwood(:,1)';
A_deterministic_kirkwood_second_moment=X_deterministic_kirkwood(:,2)';


%Make figures for each of the different moment closure approximations
f1=figure;
f2=figure;
f3=figure;
f4=figure;
f5=figure;

%Plot the first 10 trajectories for A
if M>=10
    %     Change to figure f1
    figure(f1)
    plot(time_vec,rec_vector_A(:,1:10));
    figure(f2)
    plot(time_vec,rec_vector_A(:,1:10));
    figure(f3)
    plot(time_vec,rec_vector_A(:,1:10));
    figure(f4)
    plot(time_vec,rec_vector_A(:,1:10));
    figure(f5)
    plot(time_vec,rec_vector_A(:,1:10));
else
    %     Change to figure f1
    figure(f1)
    plot(time_vec,rec_vector_A);
    figure(f2)
    plot(time_vec,rec_vector_A);
    figure(f3)
    plot(time_vec,rec_vector_A);
    figure(f4)
    plot(time_vec,rec_vector_A);
    figure(f5)
    plot(time_vec,rec_vector_A);
end

% Change to figure f1
figure(f1)
hold on
%Also plot on the mean over 1000 realisations
[hA1]=plot(time_vec,mean_A,'r','linewidth',5);

%Also plot the mean-field approximation to the mean
[hA2]=plot(t_deterministic_mean_field,A_deterministic_mean_field,'k--','linewidth',5);

%Change the axes to make sure they are all consistent
axis([0 T_final 0 A_init])

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of A particles')

legend([hA1,hA2],['Mean of M=',num2str(M),' simulations'],'Solution of mean-field model');

% Change to figure f2
figure(f2)
hold on
%Also plot on the mean over M realisations
[hA1]=plot(time_vec,mean_A,'r','linewidth',5);

%Also plot the pairwise approximation to the mean
[hA2]=plot(t_deterministic_pairwise,A_deterministic_pairwise,'k--','linewidth',5);

[hA3]=plot(time_vec,mean_A+2*sqrt(var_A),'b','linewidth',5);
plot(time_vec,mean_A-2*sqrt(var_A),'b','linewidth',5)

%Calculate the pairwise standard deviation
SD_pairwise=sqrt(A_deterministic_pairwise_second_moment-A_deterministic_pairwise.^2);

%Also plot the pairwise approxmation to +/- 2 SDs
[hA4]=plot(t_deterministic_pairwise,A_deterministic_pairwise+2*SD_pairwise,'g--','linewidth',5);
plot(t_deterministic_pairwise,A_deterministic_pairwise-2*SD_pairwise,'g--','linewidth',5)

%Change the axes to make sure they are all consistent
axis([0 T_final 0 A_init])


% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of A particles')

legend([hA1,hA3,hA2,hA4],['Mean of M=',num2str(M),' simulations'],'Mean of simulations +/- 2 SD'...
    ,'Mean of pairwise approximation','Mean of pairwise approximation +/- 2 SD');

% Change to figure f3
figure(f3)
hold on
%Also plot on the mean over M realisations
[hA1]=plot(time_vec,mean_A,'r','linewidth',5);

%Also plot the pairwise approximation to the mean
[hA2]=plot(t_deterministic_normal,A_deterministic_normal,'k--','linewidth',5);

[hA3]=plot(time_vec,mean_A+2*sqrt(var_A),'b','linewidth',5);
plot(time_vec,mean_A-2*sqrt(var_A),'b','linewidth',5)

%Calculate the pairwise standard deviation
SD_normal=sqrt(A_deterministic_normal_second_moment-A_deterministic_normal.^2);

%Also plot the pairwise approxmation to +/- 2 SDs
[hA4]=plot(t_deterministic_normal,A_deterministic_normal+2*SD_normal,'g--','linewidth',5);
plot(t_deterministic_normal,A_deterministic_normal-2*SD_normal,'g--','linewidth',5)

%Change the axes to make sure they are all consistent
axis([0 T_final 0 A_init])


% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of A particles')

legend([hA1,hA3,hA2,hA4],['Mean of M=',num2str(M),' simulations'],'Mean of simulations +/- 2 SD'...
    ,'Mean of normal approximation','Mean of normal approximation +/- 2 SD');


% Change to figure f4
figure(f4)
hold on
%Also plot on the mean over M realisations
[hA1]=plot(time_vec,mean_A,'r','linewidth',5);

%Also plot the pairwise approximation to the mean
[hA2]=plot(t_deterministic_fourth_cumulant,A_deterministic_fourth_cumulant,'k--','linewidth',5);

[hA3]=plot(time_vec,mean_A+2*sqrt(var_A),'b','linewidth',5);
plot(time_vec,mean_A-2*sqrt(var_A),'b','linewidth',5)

%Calculate the pairwise standard deviation
SD_fourth_cumulant=sqrt(A_deterministic_fourth_cumulant_second_moment-A_deterministic_fourth_cumulant.^2);

%Also plot the pairwise approxmation to +/- 2 SDs
[hA4]=plot(t_deterministic_fourth_cumulant,A_deterministic_fourth_cumulant+2*SD_fourth_cumulant,'g--','linewidth',5);
plot(t_deterministic_fourth_cumulant,A_deterministic_fourth_cumulant-2*SD_fourth_cumulant,'g--','linewidth',5)

%Change the axes to make sure they are all consistent
axis([0 T_final 0 A_init])


% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of A particles')

legend([hA1,hA3,hA2,hA4],['Mean of M=',num2str(M),' simulations'],'Mean of simulations +/- 2 SD'...
    ,'Mean of 4th cumulant approximation','Mean of 4th cumulant approximation +/- 2 SD');

% Change to figure f5
figure(f5)
hold on
%Also plot on the mean over M realisations
[hA1]=plot(time_vec,mean_A,'r','linewidth',5);

%Also plot the pairwise approximation to the mean
[hA2]=plot(t_deterministic_kirkwood,A_deterministic_kirkwood,'k--','linewidth',5);

[hA3]=plot(time_vec,mean_A+2*sqrt(var_A),'b','linewidth',5);
plot(time_vec,mean_A-2*sqrt(var_A),'b','linewidth',5)

%Calculate the pairwise standard deviation
SD_kirkwood=sqrt(A_deterministic_kirkwood_second_moment-A_deterministic_kirkwood.^2);

%Also plot the pairwise approxmation to +/- 2 SDs
[hA4]=plot(t_deterministic_kirkwood,A_deterministic_kirkwood+2*SD_kirkwood,'g--','linewidth',5);
plot(t_deterministic_kirkwood,A_deterministic_kirkwood-2*SD_kirkwood,'g--','linewidth',5)

%Change the axes to make sure they are all consistent
axis([0 T_final 0 A_init])

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of A particles')


legend([hA1,hA3,hA2,hA4],['Mean of M=',num2str(M),' simulations'],'Mean of simulations +/- 2 SD'...
    ,'Mean of Kirkwood approximation','Mean of Kirkwood approximation +/- 2 SD');


%     Change to figure f1
figure(f1)
exportfig(gcf,...
    ['dimerisation_mean_field_moment_closure.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',17);
%             'LineWidth',6);
%Save as a .fig as well
saveas(gcf,['dimerisation_mean_field_moment_closure.fig'],'fig');

%     Change to figure f2
figure(f2)
exportfig(gcf,...
    ['dimerisation_pairwise_moment_closure.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',17);
%             'LineWidth',6);
%Save as a .fig as well
saveas(gcf,['dimerisation_pairwise_moment_closure.fig'],'fig');

%     Change to figure f3
figure(f3)
exportfig(gcf,...
    ['dimerisation_normal_moment_closure.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',17);
%             'LineWidth',6);
%Save as a .fig as well
saveas(gcf,['dimerisation_normal_moment_closure.fig'],'fig');

% Change to figure f4
figure(f4)
exportfig(gcf,...
    ['dimerisation_fourth_cumulant_moment_closure.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',17);
%             'LineWidth',6);
%Save as a .fig as well
saveas(gcf,['dimerisation_fourth_cumulant_moment_closure.fig'],'fig');

% Change to figure f5
figure(f5)
exportfig(gcf,...
    ['dimerisation_kirkwood_moment_closure.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',17);
%             'LineWidth',6);
%Save as a .fig as well
saveas(gcf,['dimerisation_kirkwood_moment_closure.fig'],'fig');