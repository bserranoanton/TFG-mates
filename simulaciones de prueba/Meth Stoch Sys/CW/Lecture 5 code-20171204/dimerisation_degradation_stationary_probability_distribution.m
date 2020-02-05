%This code is desgined to simulate the stochatic dimerisation degradation
%and production reactions (4.1) and to calculate the stationary probability distribution.
% This will reproduce figure 4.1 (b)
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
T_final=500;

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
                %Implement the production reaction
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

%Define the time after which we consider the process to be stationary
t_stat=30;

%Disregard the first t_stat seconds worth of data
rec_vector_A_cut=rec_vector_A(t_stat/rec_step:end,:);

%Set the manximum value of X
max_data_A=max(rec_vector_A_cut(:));

%Decide the range into which the data will fall
X_A=[0:1:max_data_A];

%extract histogram data from each of the repeats after 40 seconds when the
%process should have roughly reached steady state
[N,X_A]=hist(rec_vector_A_cut(:),X_A);

%Normalised the data to make it a PDF
N=N/sum(N);

% Plot the histogram of the stationary frequencies
[h1]=bar(X_A,N);

%correct the axes
axis([0 max_data_A+0.5 0 ceil(max(N)*100)/100])

% First calculate the normalisation constant C
C=1/(sqrt(2)*besseli(1,2*sqrt(2*k2nu/k1onu)));

% Define a vector for the stationary probability distribution a priori
phi=zeros(1,max_data_A+1);

% Then for all the values required calculate the stationary probability distribution
for i=0:max_data_A
    phi(i+1)=(C/factorial(i))*((k2nu/k1onu)^(i/2))*besseli(i-1,2*sqrt(k2nu/k1onu));
end

hold on
%Also plot the analytically determined mean
[h2]=plot(X_A,phi,'r','linewidth',5);

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('number of particles')
ylabel('stationary distribution')

legend([h1,h2],['Long-run SSA'],'PGF analysis');

exportfig(gcf,...
            ['dimerisation_degradation_SPD.eps'],...
            'Format','eps2',...
            'Width','20',...
            'Color','cmyk',...
            'Resolution',300,...
            'FontMode','fixed',...
            'FontSize',21);
%             'LineWidth',6);
        %Save as a .fig as well
        saveas(gcf,['dimerisation_degradation_SPD.fig'],'fig');