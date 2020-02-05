%This code is desgined to simulate the stochatic resonance example of
%lecture 7  defined by reactions (7.1)
%By Kit Yates
%Created 05/011/15
%Last Modified 05/11/15

clear all
close all

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

%Define speed_factor which will tell us ho much to increase the parameter
%values by to simulate in minutes rather than seconds
speed_factor=60;

%Define the rate constant of degradation
k1onu2=4*10^(-5)*speed_factor; %(sec^{-1})
k2nu=50*speed_factor; %(sec^{-1})
k3=10*speed_factor; %(sec^{-1})

if PARAMS_1
    k4nu=25*speed_factor; %(sec^{-1})
end

if PARAMS_2
    k4nu=100*speed_factor; %(sec^{-1})
end

% Define the number of repeats we will do
M=1;

% Define the final time we will simulate to
T_final=80;

%Define the initial number of particles
A_init=10;
B_init=10;


%Define how long the recording vectors should be
num_rec_steps=10^6;

%Initialise the index which will tell us where to write the current values
rec_ind=1;

%Instantiate a vector which will hold the time varying values of A
rec_vector_A=-ones(1,num_rec_steps);
rec_vector_B=-ones(1,num_rec_steps);

%Write the initial condisiotn to these vectors
rec_vector_A(rec_ind)=A_init;
rec_vector_B(rec_ind)=B_init;

% Initialise a vector which will hold the times when reactions occur
time_vec=-ones(1,num_rec_steps);

%Write the initial time to this vector
time_vec(rec_ind)=0;


%Define a vector which will hold the propensity functions
a=zeros(4,1);

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
    
    %Increase the recording index
    rec_ind=rec_ind+1;
    
    %Define the propensity functions
    a(1) =  k1onu2 * A* (A-1)*B;
    a(2) =  k2nu;
    a(3) =  k3 * A;
    a(4) =  k4;
    
    %Calculate the cumulative sum of a
     cumsuma=cumsum(a);
    
    %Calculate the propensity functions
    a0 = cumsuma(end);
    
    %Determine the time for the next reaction
    tau =(1/a0)*log(1/rand);
    
    %Update the time
    t = t + tau;
    
    %Write the current value of A to A_old
    A_old=A;
    B_old=B;
    
    %multiply the random number by the sum of the propensities and
    %store so we can reuse
    ra0=rand*a0;
    
    %if we have not run out of molecules
    if tau<inf;
        %Determine which reaction happened
        if ra0 < cumsuma(1)
            %Implement the thrd order production reaction
            A = A_old  +1;
            B = B_old -1;
        elseif ra0 < cumsuma(2)
            %Implement the zeroth order production reaction
            A = A_old +1;
        elseif ra0 < cumsuma(3)    
            % Implement a first order degradation reaction
            A = A_old -1;
        else
            % Implement the zeroth order poduction reaction for B
            B = B_old +1;
        
        end
    end
    
    %Rcord the time and the numbers of molecules
    time_vec(rec_ind)=t;
    rec_vector_A(rec_ind)=A;
    rec_vector_B(rec_ind)=B;
    
end


%Get rid of the residual -1s in the recording vector
time_vec(time_vec<0)=[];
rec_vector_A(rec_vector_A<0)=[];
rec_vector_B(rec_vector_B<0)=[];

%Solve the ODEs for the deterministic system
[t_deterministic,X_deterministic]=ode15s(@(t,X)RHS(t,X,k1onu2,k2nu,k3,k4nu),[0,T_final],[A_init,B_init]);

%Write A and B from X
A_deterministic=X_deterministic(:,1)';
B_deterministic=X_deterministic(:,2)';


if PARAMS_1
    figure
    %plot the stochastic trajectory for A
    [hA1]=semilogy(time_vec,rec_vector_A,'r');
    hold on
    %plot the determinisitc trajectory fo A
    [hA2]=semilogy(t_deterministic,A_deterministic,'k--','linewidth',5);
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('number of A particles')
    
    %Find the original axis
    original_axis=axis;
    
    %Make sure the axes only go up to th final time
    axis([0 80 original_axis(3) original_axis(4)])
    
    
    legend([hA1,hA2],'Stochastic trajectory','Deterministic trajectory');
    
    %
    exportfig(gcf,...
        ['resonance_stochastic_deterministic_A_log_scale.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['resonance_stochastic_deterministic_A_log_scale.fig'],'fig');
    
    %   Make a new figure
    figure
    %plot the stochastic trajectory for A
    [hA]=plot(time_vec,rec_vector_A,'r','linewidth',5);
    hold on
    %don't plot all the values of B so we can see A.
    %     Decide how often to plot B
    B_refinement_factor=5000;
    %plot the determinisitc trajectory for A
    [hB]=plot(time_vec(1:B_refinement_factor:end),rec_vector_B(1:B_refinement_factor:end),'g--','linewidth',5);
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('number of particles')
    
    %Find the original axis
    original_axis=axis;
    
    %Make sure the axes only go up to th final time
    axis([0 80 original_axis(3) original_axis(4)])
    
    legend([hA,hB],'A(t)','B(t)');
    
    %
    exportfig(gcf,...
        ['resonance_stochastic_A_and_B_linear_scale.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['resonance_stochastic_A_and_B_linear_scale.fig'],'fig');
    
    figure
    %Plot the stochastic vs detemrinistic phase planes and nullclines
    %plot the stochastic trajectory
    [hAB_stoch]=semilogx(rec_vector_A,rec_vector_B,'r');
    hold on
    original_axis=axis;
    %plot the nullcline for da/dt=0
    a_for_null=1:1:10^4;
    b_for_null_a=(k3*a_for_null-k2nu)./(a_for_null.^2*k1onu2);
    [h_null_for_a]=semilogx(a_for_null,b_for_null_a,'g','linewidth',5);
    %plot the nullcline for db/dt=0
    b_for_null_b=k4nu./(a_for_null.^2*k1onu2);
    [h_null_for_b]=semilogx(a_for_null,b_for_null_b,'g','linewidth',5);
    %Plot the deterministic trajectory
    [hAB_det]=semilogx(A_deterministic,B_deterministic,'k--','linewidth',3);
    %Finally Plot the  intercept point
    %Calculate the value of A
    a_intercept=(k2nu+k4nu)/k3;
    %calculate the value of B
    b_intercept=(k3*a_intercept-k2nu)./(a_intercept.^2*k1onu2);
    %Plot the point as a blue dot
    semilogx(a_intercept,b_intercept,'bo','markersize',5,'linewidth',5)
    %Update the axis
    axis([original_axis(1) original_axis(2) original_axis(3) max(b_for_null_a)])
    %Add labels to axes
    xlabel('number of A particles')
    ylabel('number of B particles')
        %
    exportfig(gcf,...
        ['resonance_phase_plane_intermediate_regime.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['resonance_phase_plane_intermediate_regime.fig'],'fig');
end

if PARAMS_2
    figure
    %plot the stochastic trajectory for A
    [hA1]=semilogy(time_vec(time_vec<20),rec_vector_A(time_vec<20),'r');
    hold on
    %plot the determinisitc trajectory fo A
    [hA2]=semilogy(t_deterministic(t_deterministic<20),A_deterministic(t_deterministic<20),'k--','linewidth',1);
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('number of A particles')
    
    %Find the original axis
    original_axis=axis;
    
    %Make sure the axes only go up to th final time
    axis([0 20 original_axis(3) original_axis(4)])
    
    
    legend([hA1,hA2],'Stochastic trajectory','Deterministic trajectory');
    
    %
    exportfig(gcf,...
        ['resonance_stochastic_deterministic_oscillatory_regime_A.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['resonance_stochastic_deterministic_oscillatory_regime_A.fig'],'fig');
    
    %   Make a new figure
    figure
    %plot the stochastic trajectory for A
    [hB1]=plot(time_vec(time_vec<20),rec_vector_B(time_vec<20),'r','linewidth',1);
    hold on

    %plot the determinisitc trajectory for A
    [hB2]=plot(t_deterministic(t_deterministic<20),B_deterministic(t_deterministic<20),'k--','linewidth',1);
    %Set the x and y labels
    xlabel('time (min)')
    ylabel('number of particles')
    
    %Find the original axis
    original_axis=axis;
    
    %Make sure the axes only go up to th final time
    axis([0 20 original_axis(3) original_axis(4)])
    
    legend([hB1,hB2],'Stochastic trajectory','Deterministic trajectory');
    
    %
    exportfig(gcf,...
        ['resonance_stochastic_deterministic_oscillatory_regime_B.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['resonance_stochastic_deterministic_oscillatory_regime_B.fig'],'fig');
    
    figure
    %Plot the stochastic vs detemrinistic phase planes and nullclines
    %plot the stochastic trajectory
    [hAB_stoch]=semilogx(rec_vector_A,rec_vector_B,'r');
    hold on
    original_axis=axis;
    %plot the nullcline for da/dt=0
    a_for_null=1:1:10^4;
    b_for_null_a=(k3*a_for_null-k2nu)./(a_for_null.^2*k1onu2);
    [h_null_for_a]=semilogx(a_for_null,b_for_null_a,'g','linewidth',5);
    %plot the nullcline for db/dt=0
    b_for_null_b=k4nu./(a_for_null.^2*k1onu2);
    [h_null_for_b]=semilogx(a_for_null,b_for_null_b,'g','linewidth',5);
    %Plot the deterministic trajectory
    [hAB_det]=semilogx(A_deterministic,B_deterministic,'k--','linewidth',3);
    %Finally Plot the  intercept point
    %Calculate the value of A
    a_intercept=(k2nu+k4nu)/k3;
    %calculate the value of B
    b_intercept=(k3*a_intercept-k2nu)./(a_intercept.^2*k1onu2);
    %Plot the point as a blue dot
    semilogx(a_intercept,b_intercept,'bo','markersize',5,'linewidth',5)
    %Update the axis
    axis([original_axis(1) max(A_deterministic) original_axis(3) max(B_deterministic)])
    %Add labels to axes
    xlabel('number of A particles')
    ylabel('number of B particles')
        %
    exportfig(gcf,...
        ['resonance_phase_plane_oscillatory_regime.eps'],...
        'Format','eps2',...
        'Width','20',...
        'Color','cmyk',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',21);
    %             'LineWidth',6);
    %Save as a .fig as well
    saveas(gcf,['resonance_phase_plane_oscillatory_regime.fig'],'fig');
end