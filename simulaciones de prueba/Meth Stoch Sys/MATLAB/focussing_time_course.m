%This code is desgined to simulate the stochastic focussing reactions in
%Lecture 5
%By Kit Yates
%Created 28/08/15
%Last Modified 26/10/16

clear all
close all

%Set up some options which allow us to change things about and plot the
%solutions

%NORMAL tells us to use the normal parameter values
NORMAL = 1;
%ALTERED tells us to use the altered parameter values
ALTERED = 0;

%AFIXED does not allow the values of A to change during the stochastic
%simulation
AFIXED=0; %a fixed

%Check that at least one of NORMAL and ALTERED is one
if (NORMAL+ALTERED)<1
    error('You have not selected a set of parameters')
end

%Check that only one of them is 1
if (NORMAL+ALTERED)>1
    error('You hve chosen too many sets of parameters')
end

%Check that if AFIXED is on then we are using the normal parameter values

if AFIXED
    if ~NORMAL
        error('You are fixing A with the wrong set of parameters')
    end
end

%
%Decide the factor by which to speed up the simulation
speed_factor=60;

%Define the rate constant of the reactions
if NORMAL
    k1nu=10^(2)*speed_factor; %(min^{-1})
    k2=10^(3)*speed_factor; %(sec^{-1})
    k3=10^(-2)*speed_factor; %(sec^{-1})
    k4onu=9900*speed_factor; %(sec^{-1})
    k5nu=10^3*speed_factor;
    k6=10^2*speed_factor;
end

if ALTERED
    k1nu=10^(2)*speed_factor; %(sec^{-1})
    k2=0.1*speed_factor; %(sec^{-1})
    k3=10^(-2)*speed_factor; %(sec^{-1})
    k4onu=0.99*speed_factor; %(sec^{-1})
    k5nu=10^3*speed_factor;
    k6=10^2*speed_factor;
end

% Define the number of repeats we will do
M=100;

% Define the final time we will simulate to
T_final = 35;

%Define the time at which the parameters change
T_change = 10;

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
rec_vector_C=zeros(num_rec_steps+1,M);

%Define the initial number of particles
A_init = 10;
B_init = 100;
C_init = 0;

%Define the initial number of molecules for the recording vector
rec_vector_A(1,:) = A_init * ones(1,M);
rec_vector_B(1,:) = B_init * ones(1,M);
rec_vector_C(1,:) = C_init * ones(1,M);

%Define a vector which will hold the propensity functions
a = zeros(6,1);

%Run through a for loop for each of the repeats
for i=1:M
    
    %Define the initial time to be zero
    t=0;
    
    %Initialise the times to help recording
    t_after=0;
    
    %   initialise the number of particles for this repeat
    A = A_init;
    B = B_init;
    C = C_init;
    
    %for each repeat initialise the variable which tells us whether we have
    %not changed the parameters yet
    not_changed=1;
    
    %Reset the parameter k5nu every time we do anothe repeat
    k5nu=10^3*speed_factor;
    
    %Run through a for loop for each of the time-steps
    %Enter the while loop
    while t<T_final
        
        %Define the propensity functions
        a(1) = k1nu;
        a(2) = k2 * C;
        a(3) = k3 * B;
        a(4) = k4onu * A * C;
        %if A is fixed then set the parameters for these reactions to zero
        if AFIXED
            a(5)=0;
            a(6)=0;
        else
            a(5)=k5nu;
            a(6) = k6*A;
        end
        
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
        C_old=C;
        
        %multiply the random number by the sum of the propensities and
        %store so we can reuse
        ra0=rand*a0;
        
        %if we have not run out of molecules
        if tau<inf
            %Determine which reaction happened
            if ra0 < cumsuma(1) 
                %Implement production of C
                C = C_old+1;
            elseif ra0 < cumsuma(2)
                %Implement the conversion of C to B
                C = C_old -1;
                B = B_old +1;
            elseif ra0 < cumsuma(3)
                % Implement a degradation of B
                B = B_old -1;
            elseif ra0 < cumsuma(4)
                % Implement degradation of C catalyse by A
                C = C_old -1;
                
            elseif ra0 < cumsuma(5)
                % Implement production of A
                A = A_old +1;
            else
                %Implement degradation of A
                A = A_old -1;
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
            rec_vector_C(ind_before+1:ind_after+1,i)=(C_old)*ones(steps_to_write,1);
        end
        
        %If we have not yet changed the parameter
        if not_changed
            %Check to see if we have reached the appropriate time
            if t > T_change
                %                 If we have then changethe parameter
                k5nu=5*10^2*speed_factor; %to convert into min
                %If A is fixed
                if AFIXED
                   %change the value of A so that it is again at steady
                   %state
                   A = k5nu/k6;
                end
                
                %and change the variable which tells us whether or not we
                %have changed
                not_changed=0;
            end
        end
        
        
    end
    
end

%Calculate the mean over all 100 realisations
mean_A=mean(rec_vector_A,2);
mean_B=mean(rec_vector_B,2);
mean_C=mean(rec_vector_C,2);

%Reset k5nu to the correct value
k5nu=10^3*speed_factor;

%Solve the ODEs for the deterministic system until we change the parameter
%values


%Change the value of k5nu


%Solve the ODEs for the deterministic system after we change the parameter
%values (if A is fixed we have to ensure )
if AFIXED
    
    
else

    
end

%Write A and B from X
A_deterministic=[X_deterministic_first(:,1)',X_deterministic_second(:,1)'];
B_deterministic=[X_deterministic_first(:,2)',X_deterministic_second(:,2)'];
C_deterministic=[X_deterministic_first(:,3)',X_deterministic_second(:,3)'];

%Concatonate the time vectors
t_deterministic=[t_deterministic_first;t_deterministic_second];

%If we are considering the normal parameters with the variable value of A
if NORMAL && ~AFIXED
    
    %Plot the first 5 trajectories for A and B
    f1=figure;
    f2=figure;
    f3=figure;
    
    if M>=10
        %     Change to figure f1
        figure(f1)
        plot(time_vec,rec_vector_A(:,1:10));
        
        %     Change to figure f2
        figure(f2)
        plot(time_vec,rec_vector_B(:,1:10));
        
                %     Change to figure f3
        figure(f3)
        plot(time_vec,rec_vector_C(:,1:10));
    else
        
        %     Change to figure f1
        figure(f1)
        plot(time_vec,rec_vector_A);
        
        %     Change to figure f2
        figure(f2)
        plot(time_vec,rec_vector_B);
                %     Change to figure f3
        figure(f3)
        plot(time_vec,rec_vector_C);
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
    xlabel('time (min)')
    ylabel('number of A particles')
    
    legend([hA1,hA2],['Mean of M=',num2str(M),' simulations'],'Solution of ODEs');
    
    
    %     Change to figure f2
    figure(f2)
    hold on
    %Also plot on the mean over 100 realisations
    [hB1]=plot(time_vec,mean_B,'r','linewidth',5);
    
    %Also plot the analytically determined mean
    [hB2]=plot(t_deterministic,B_deterministic,'k--','linewidth',5);
    
    % %Set the title
    % title(['t= ', num2str(t*rec_step_increment)])
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('number of B particles')
    
    legend([hB1,hB2],['Mean of M=',num2str(M),' simulations'],'Solution of ODEs');
    
        %     Change to figure f3
    figure(f3)
    hold on
    %Also plot on the mean over 100 realisations
    [hC1]=plot(time_vec,mean_C,'r','linewidth',5);
    
    %Also plot the analytically determined mean
    [hC2]=plot(t_deterministic,C_deterministic,'k--','linewidth',5);
    
    % %Set the title
    % title(['t= ', num2str(t*rec_step_increment)])
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('number of C particles')
    
    legend([hC1,hC2],['Mean of M=',num2str(M),' simulations'],'Solution of ODEs');
    
    
    %     Change to figure f1
    figure(f1)
    exportfig(gcf,...
        ['focussing_ODE_and_paths_for_A.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['focussing_ODE_and_paths_for_A.fig'],'fig');
    
    %     Change to figure f2
    figure(f2)
    exportfig(gcf,...
        ['focussing_ODE_and_paths_for_B.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['focussing_ODE_and_paths_for_B.fig'],'fig');
        
    %     Change to figure f2
    figure(f3)
    exportfig(gcf,...
        ['focussing_ODE_and_paths_C.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['focussing_ODE_and_paths_for_C.fig'],'fig');
    
end

if NORMAL && AFIXED
   %Plot the first 5 trajectories for A and B
    f1=figure;
    f2=figure;
    
    if M>=10
        %     Change to figure f1
        figure(f1)
        plot(time_vec,rec_vector_A(:,1:10));
        
        %     Change to figure f2
        figure(f2)
        plot(time_vec,rec_vector_B(:,1:10));
    else
        
        %     Change to figure f1
        figure(f1)
        plot(time_vec,rec_vector_A);
        
        %     Change to figure f2
        figure(f2)
        plot(time_vec,rec_vector_B);
    end
    
    %     Change to figure f1
    figure(f1)
    hold on
    %Also plot on the mean over 100 realisations
    [hA1]=plot(time_vec,mean_A,'r','linewidth',5);
    
    %Also plot the analytically determined mean
    [hA2]=plot(t_deterministic,A_deterministic,'k--','linewidth',5);
    
    %Change this axes
    axis([0 35 0 25])
    
    % %Set the title
    % title(['t= ', num2str(t*rec_step_increment)])
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('number of A particles')
    
    legend([hA1,hA2],['Mean of M=',num2str(M),' simulations'],'Solution of ODEs');
    
    
    %     Change to figure f2
    figure(f2)
    hold on
    %Also plot on the mean over 100 realisations
    [hB1]=plot(time_vec,mean_B,'r','linewidth',5);
    
    %Also plot the analytically determined mean
    [hB2]=plot(t_deterministic,B_deterministic,'k--','linewidth',5);
    
    % %Set the title
    % title(['t= ', num2str(t*rec_step_increment)])
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('number of B particles')
    
    legend([hB1,hB2],['Mean of M=',num2str(M),' simulations'],'Solution of ODEs','location','SouthEast');
    
    %     Change to figure f1
    figure(f1)
    exportfig(gcf,...
        ['focussing_A_fixed.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['focussing_A_fixed.fig'],'fig');
    
    %     Change to figure f2
    figure(f2)
    exportfig(gcf,...
        ['Focussing_ODE_and_Paths_B_A_fixed.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['Focussing_ODE_and_Paths_B_A_fixed.fig'],'fig');
         
end

if ALTERED
    
    %Plot the first 5 trajectories for A and B
    f1=figure;
    f2=figure;
    f3=figure;
    
    if M>=10
        %     Change to figure f1
        figure(f1)
        plot(time_vec,rec_vector_A(:,1:10));
        
        %     Change to figure f2
        figure(f2)
        plot(time_vec,rec_vector_B(:,1:10));
        
                %     Change to figure f3
        figure(f3)
        plot(time_vec,rec_vector_C(:,1:10));
    else
        
        %     Change to figure f1
        figure(f1)
        plot(time_vec,rec_vector_A);
        
        %     Change to figure f2
        figure(f2)
        plot(time_vec,rec_vector_B);
                %     Change to figure f3
        figure(f3)
        plot(time_vec,rec_vector_C);
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
    xlabel('time (min)')
    ylabel('number of A particles')
    
    legend([hA1,hA2],['Mean of M=',num2str(M),' simulations'],'Solution of ODEs');
    
    
    %     Change to figure f2
    figure(f2)
    hold on
    %Also plot on the mean over 100 realisations
    [hB1]=plot(time_vec,mean_B,'r','linewidth',5);
    
    %Also plot the analytically determined mean
    [hB2]=plot(t_deterministic,B_deterministic,'k--','linewidth',5);
    
    % %Set the title
    % title(['t= ', num2str(t*rec_step_increment)])
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('number of B particles')
    
    legend([hB1,hB2],['Mean of M=',num2str(M),' simulations'],'Solution of ODEs');
    
        %     Change to figure f3
    figure(f3)
    hold on
    %Also plot on the mean over 100 realisations
    [hC1]=plot(time_vec,mean_C,'r','linewidth',5);
    
    %Also plot the analytically determined mean
    [hC2]=plot(t_deterministic,C_deterministic,'k--','linewidth',5);
    
    % %Set the title
    % title(['t= ', num2str(t*rec_step_increment)])
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('number of C particles')
    
    legend([hC1,hC2],['Mean of M=',num2str(M),' simulations'],'Solution of ODEs');
    
    
    %     Change to figure f1
    figure(f1)
    exportfig(gcf,...
        ['focussing_ODEs_and_paths_for_A_new_params.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['focussing_ODEs_and_paths_for_A_new_params.fig'],'fig');
    
    %     Change to figure f2
    figure(f2)
    exportfig(gcf,...
        ['focussing_ODEs_and_paths_for_B_new_params.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['focussing_ODEs_and_paths_for_B_new_params.fig'],'fig');
        
    %     Change to figure f2
    figure(f3)
    exportfig(gcf,...
        ['focussing_ODEs_and_paths_for_C_new_params.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['focussing_ODEs_and_paths_for_C_new_params.fig'],'fig'); 
    
end