%This code is desgined to simulate the Schlogl reaction system (6.1).
% This will reproduce figures 6.1 and 6.2 (a)
%By Kit Yates
%Created 19/10/15
%Last Modified 19/10/15

clear all
close all

%Define some options
%This option solves the PDE with mutliple ICs and plots the corresponding
%figure
MULTIPLE_ICs=0;

%This initial condition plots the stochast and deterministic system
%starting from a zero IC
ZERO_IC=1;

%Check that at least one of NORMAL and ALTERED is one
if (MULTIPLE_ICs+ZERO_IC)<1
    error('You have not selected a set of initial conditions')
end

%Check that only one of them is 1
if (MULTIPLE_ICs+ZERO_IC)>1
    error('You hve chosen too many sets of initial conditions')
end

%If we have decided to initialise with multiple ICs then decide what they
%are
IC_list=[500,300,200,0];

%Find the number of initial conditions
num_ICs=length(IC_list);

%Define the rate constants
k1onu2=2.5*10^(-4); %(min^{-1})
k2onu=0.18; %(min^{-1})
k3=37.5; %(min^{-1})
k4nu=2200; %(min^{-1})

% Define the final time we will simulate to
T_final=100;

%Only run the stochastic simulation in the case when we are plotting a
%single deterministic trajectory to compare against
if ZERO_IC
    % Define the number of repeats we will do
    M=1;
    
    
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
    a=zeros(4,1);
    
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
            a(1)=A*(A-1)*(A-2)*k1onu2;
            a(2)=A*(A-1)*k2onu;
            a(3)=A*k3;
            a(4)=k4nu;
            
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
            
            %multiply the random number by the sum of the propensities and
            %store so we can reuse
            ra0=rand*a0;
            
            %if we have not run out of molecules
            if tau<inf;
                %Determine which reaction happened
                if ra0<cumsuma(1)
                    %Implement the trimerisation degradation reaction
                    A=A_old-1;
                elseif ra0<cumsuma(2)
                    % Implement a bimolecular production reaction
                    A=A_old+1;
                elseif ra0<cumsuma(3)
                    %Implement the degradation reaction
                    A=A_old-1;
                else
                    %Implement the production reaction
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
    if M==1
        mean_A=rec_vector_A;
    else
        mean_A=mean(rec_vector_A,2);
    end
    
    
end
if MULTIPLE_ICs
        %Solve the ODEs for the deterministic system
        [t_deterministic_1,A_deterministic_1]=ode15s(@(t,X)RHS_bistable(t,X,k1onu2,k2onu,k3,k4nu),[0,1],IC_list(1)); 
                %Solve the ODEs for the deterministic system
        [t_deterministic_2,A_deterministic_2]=ode15s(@(t,X)RHS_bistable(t,X,k1onu2,k2onu,k3,k4nu),[0,1],IC_list(2));
                %Solve the ODEs for the deterministic system
        [t_deterministic_3,A_deterministic_3]=ode15s(@(t,X)RHS_bistable(t,X,k1onu2,k2onu,k3,k4nu),[0,1],IC_list(3));
                %Solve the ODEs for the deterministic system
        [t_deterministic_4,A_deterministic_4]=ode15s(@(t,X)RHS_bistable(t,X,k1onu2,k2onu,k3,k4nu),[0,1],IC_list(4));
end

% Solv the ODE for the single zero IC
if ZERO_IC
    %Solve the ODEs
    [t_deterministic,A_deterministic]=ode15s(@(t,X)RHS_bistable(t,X,k1onu2,k2onu,k3,k4nu),[0,T_final],0);
end

if MULTIPLE_ICs
    %plot the positions of the two stable teady states
    plot([0:0.01:1],400,'k--','linewidth',3)
    hold on
    plot([0:0.01:1],100,'k--','linewidth',3)
    %Plot the deterministic solution
    [h_500]=plot(t_deterministic_1,A_deterministic_1,'b','linewidth',5);
    %Plot the deterministic solution
    [h_300]=plot(t_deterministic_2,A_deterministic_2,'g','linewidth',5);
    %Plot the deterministic solution
    [h_200]=plot(t_deterministic_3,A_deterministic_3,'k','linewidth',5);
    %Plot the deterministic solution
    [h_0]=plot(t_deterministic_4,A_deterministic_4,'r','linewidth',5);
    
    
        %lable the axes
    xlabel('time (min)')
    ylabel('number of particles')
    
        %add a legend
    legend([h_500,h_300,h_200,h_0],['A(0)=',num2str(IC_list(1))],['A(0)=',num2str(IC_list(2))],...
        ['A(0)=',num2str(IC_list(3))],['A(0)=',num2str(IC_list(4))],'location','East');
    
%     Save the figure
    exportfig(gcf,...
        ['schlogl_deterministic_evolutions.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['schlogl_deterministic_evolutions.fig'],'fig');
    
end

if ZERO_IC
   %Define the time to crop the solutions at
    T_crop=2;
    
    time_vec_cropped=[0:rec_step:T_crop];
    mean_A_cropped=mean_A(1:length(time_vec_cropped));
    %plot the cropped stochastic trajectory
    figure
    [hA_stoch]=plot(time_vec_cropped,mean_A_cropped);
    hold on
    
    %Plot the deterministic solution on top
    [hA_det]=plot(t_deterministic(t_deterministic<T_crop),A_deterministic(t_deterministic<T_crop),'k--','linewidth',5);
    
    %lable the axes
    xlabel('time (min)')
    ylabel('number of particles')
    
    %add a legend
    legend([hA_stoch,hA_det],'stochastic','deterministic');
    
    %Save the fgure
    exportfig(gcf,...
        ['schlogl_deterministic_and_stochastic_short_time.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['schlogl_deterministic_and_stochastic_short_time.fig'],'fig');
    
    figure
    %plot the full stochastic trajectory
    [hA_stoch]=plot(time_vec,mean_A);
    hold on
    %Plot the deterministic solution on top
    [hA_det]=plot(t_deterministic,A_deterministic,'k--','linewidth',5);
    
    %lable the axes
    xlabel('time (min)')
    ylabel('number of particles')
    
    %add a legend
    legend([hA_stoch,hA_det],'stochastic','deterministic');
    
    %save the figure
    exportfig(gcf,...
        ['schlogl_deterministic_and_stochastic_long_time.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['schlogl_deterministic_and_stochastic_long_time.fig'],'fig');
end