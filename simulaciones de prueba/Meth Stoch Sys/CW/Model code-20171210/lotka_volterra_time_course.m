%This code is desgined to simulate the stochatic resonance example of
%lecture 7  defined by reactions (7.1)
%By Kit Yates
%Created 07/11/16
%Last Modified 07/11/16

clear all
close all
tic;
% Define Booleans which between which parameter values we use
PARAMS_1=1; %Determinsitic system does not oscillate
PARAMS_2=0; %Deterministic system oscillates

%%Check that we have chosen at least one set of parameter values
if (PARAMS_1+PARAMS_2)<1
    error('You have not selected a set of parameters')
end

%Check that we have not chosen more than one set of parameter values.
if (PARAMS_1+PARAMS_2)>1
    error('You hve chosen too many sets of parameters')
end

if PARAMS_1
    %Define the rate of birth of prey
    k1=3; %(sec^{-1})
    %Define the rate of conversion from predator to prey
    k2onu=0.01; %(sec^{-1})
    %Define the rate of death of predators
    k3=3; %(sec^{-1})
end


% Define the number of repeats we will do
M=1000;

%Define the number of repeats to record
M_rec=1000;

% Define the final time we will simulate to
T_final=50;

%Define the recording time interval
rec_step=0.01;

%Define a vector of the time_points
time_vec=[0:rec_step:T_final];

%Calculate the number of recording steps required
num_rec_steps=T_final/rec_step;

%Define the vector which will reocrd the number of animalsat each time
%point A is prey and B is predators
rec_vector_A=zeros(num_rec_steps+1,M_rec);
rec_vector_B=zeros(num_rec_steps+1,M_rec);

%Define the initial number of animals
A_init=(k3+1)/k2onu; %steady state is at k3/k2onu
B_init=(k1+1)/k2onu; %steady state is at k1/k2onu

%Because of the initial condition we start in the top right quadrant with
%respect to the steady state Q=1 to 4 defines which quadrant we are in.
%They are numbered anticlockwise starting with one being the upper right
%quadrant
Q=1;

% Calculate the values of the non-trivial steady state for A and B
A_fixed=k3/k2onu;
B_fixed=k1/k2onu;

%Define the initial number of molecules for the recording vector
rec_vector_A(1,:)=A_init*ones(1,M_rec);
rec_vector_B(1,:)=B_init*ones(1,M_rec);

%Define a vector which will hold the propensity functions
a=zeros(3,1);

%Initialise vecotres which will hold the number of times rey and predators
%go extict
T_prey_extinction=[];
T_predator_extinction=[];

%Define variables which will hold the sum of the periods and the sum of the
%square periods so we can calculate the sample variance easily
sum_periods=0;
sum_periods_square=0;

%Initialise a counter to say we are on the first period
period_counter=0;

%Initialise sample vairance to be very large so we enter the repeat while
%loop
var_sample_mean=inf;

%Initialise a repeat counting variable to be zero
i=0;

period_length=zeros(1,M*T_final/2);

%Run through a for loop for each of the repeats
while var_sample_mean>10^(-4) || i<M_rec;
    
    %Increase the counting index
    i=i+1;
    
    % Set the time since we last entered this top right quadrant to be zero
    % initially
    Q1_time=0;
    
    %Define the initial time to be zero
    t=0;
    
    %Initialise the times to help recording
    t_after=0;
    
    %   initialise the number of animals for this repeat
    A=A_init;
    B=B_init;
    
    %Initialise the octant of the plane we are in
    
    %Run through a for loop for each of the time-steps
    %Enter the while loop
    while t<T_final && A>0 && B>0
        
        %Define the propensity functions
        a(1)=k1*A;
        a(2)=A*B*k2onu;
        a(3)=B*k3;
        
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
                %Implement production of C
                A=A_old+1;
            elseif ra0<cumsuma(2)
                A=A_old-1;
                B=B_old+1;
            else
                B=B_old-1;
            end
        end
        %Calculate the times for recording
        t_before=t_after;
        t_after=t;
        
        
        %Determine whether the stochastic trajectory has moved quadrants
        
        if Q==1;
            %find where we are with respect to the fixed point
            if A-A_fixed<0 && B-B_fixed>0
                %Change which quadrant we are in
                Q=2;
            end
        elseif Q==2;
            %%find where we are with respect to the fixed point
            if A-A_fixed<0 && B-B_fixed<0
                %Change which quadrant we are in
                Q=3;
            end
        elseif Q==3
            %find where we are with respect to the fixed point
            if A-A_fixed>0 && B-B_fixed<0
                %Change which quadrant we are in
                Q=4;
            end
        elseif Q==4
            %find where we are with respect to the fixed point
            if A-A_fixed>0 && B-B_fixed>0
                %Change which quadrant we are in
                Q=1;
                if Q1_time~=0;
                    %sum periods
                    sum_periods=sum_periods+(t-Q1_time);
                    %Sum periods square
                    sum_periods_square=sum_periods_square+(t-Q1_time)^2;
                    %Increase the counter of the number of periods
                    period_counter=period_counter+1;
                end
                %Reset the time that we last entered the positive quadrant
                %to be the current time.
                Q1_time=t;
            end
        end
        
        
        
        if i<=M_rec
            
            %Record the extinction times if this occurs
            if A==0
                T_prey_extinction=[T_prey_extinction,t];
            end
            if B==0
                T_predator_extinction=[T_predator_extinction,t];
            end
            
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
    
    %Calculate the sample variance of the period
    var_sample_mean=(sum_periods_square/period_counter-(sum_periods/period_counter)^2)/period_counter;
    
end

%Calculate the mean over all 100 realisations
mean_A=mean(rec_vector_A,2);
mean_B=mean(rec_vector_B,2);

%Get rid of the first couple of entries of the stochastic period in each
%repeat
mean_stochastic_period=sum_periods/period_counter;

%Solve the ODEs for the deterministic system until we change the parameter
%values
[t_deterministic,X_deterministic]=ode15s(@(t,X)RHS_LV(t,X,k1,k2onu,k3),[0,T_final],[A_init,B_init]);

%Write A and B from X
A_deterministic=X_deterministic(:,1);
B_deterministic=X_deterministic(:,2);

rec_vector_A_first_repeat=rec_vector_A(:,1);
rec_vector_B_first_repeat=rec_vector_B(:,1);

% Run through a for loop for the dterministic time step and calculate the
% period
Q=1;
%Initialise a variable which will hold the number of oscillations to be
%zero
period_counter_det=0;
%initialise a vector which holds the length of the periods
period_length_det=zeros(1,T_final/2);

%Set the time of entering the first quadrant to zero
Q1_time=0;

for i=1:length(X_deterministic)
    A=A_deterministic(i);
    B=B_deterministic(i);
    t=t_deterministic(i);
    if Q==1;
            %find where we are with respect to the fixed point
            if A-A_fixed<0 && B-B_fixed>0
                %Change which quadrant we are in
                Q=2;
            end
        elseif Q==2;
            %%find where we are with respect to the fixed point
            if A-A_fixed<0 && B-B_fixed<0
                %Change which quadrant we are in
                Q=3;
            end
        elseif Q==3
            %find where we are with respect to the fixed point
            if A-A_fixed>0 && B-B_fixed<0
                %Change which quadrant we are in
                Q=4;
            end
        elseif Q==4
            %find where we are with respect to the fixed point
            if A-A_fixed>0 && B-B_fixed>0
                %Change which quadrant we are in
                Q=1;
                if Q1_time~=0;
                    period_counter_det=period_counter_det+1;
                    %Record the actual periods too
                    period_length_det(period_counter_det)=t-Q1_time;

                end
                %Reset the time that we last entered the positive quadrant
                %to be the current time.
                Q1_time=t;
            end
        end
end

%Remove any extra values of period length
period_length_det(period_length_det==0)=[];

%Calculate the mean detemrinistic period
mean_period_det=mean(period_length_det);

toc

%If we are considering the first set of parameter values
if PARAMS_1
    
    figure
    %plot the stochastic trajectory for A
    [hA1]=plot(time_vec,rec_vector_A_first_repeat,'r');
    hold on
    %plot the determinisitc trajectory fo A
    [hA2]=plot(t_deterministic,A_deterministic,'k--','linewidth',5);
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('number of prey')
    
    %Find the original axis
    original_axis=axis;
    
    %Make sure the axes only go up to th final time
    axis([0 T_final original_axis(3) original_axis(4)])
    
    
    legend([hA1,hA2],'Stochastic trajectory','Deterministic trajectory');
    
    %
    exportfig(gcf,...
        ['predator_prey_stochastic_deterministic_A_survival_regime.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['predator_prey_stochastic_deterministic_A_survival_regime.fig'],'fig');
    
    %   Make a new figure
    figure
    %plot the stochastic trajectory for A
    [hA]=plot(time_vec,rec_vector_A_first_repeat,'r','linewidth',5);
    hold on
    %plot the determinisitc trajectory for A
    [hB]=plot(time_vec,rec_vector_B_first_repeat,'g--','linewidth',5);
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('number of animals')
    
    %Find the original axis
    original_axis=axis;
    
    %Make sure the axes only go up to the final time
    axis([0 T_final original_axis(3) original_axis(4)])
    
    legend([hA,hB],'A(t)','B(t)');
    
    %
    exportfig(gcf,...
        ['predator_prey_stochastic_A_and_B_survival_regime.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['predator_prey_stochastic_A_and_B_survival_regime.fig'],'fig');
    
    figure
    %Plot the stochastic vs detemrinistic phase planes and nullclines
    %plot the stochastic trajectory
    [hAB_stoch]=plot(rec_vector_A_first_repeat,rec_vector_B_first_repeat,'r');
    hold on
    original_axis=axis;
    %Plot the deterministic trajectory
    [hAB_det]=semilogx(A_deterministic,B_deterministic,'k--','linewidth',3);
    %Finally Plot the  intercept point
    %Calculate the value of A
    %Add labels to axes
    xlabel('number of prey')
    ylabel('number of predators')
    %
    exportfig(gcf,...
        ['predator_prey_phase_plane_survival_regime.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['predator_prey_phase_plane_survival_regime.fig'],'fig');
end