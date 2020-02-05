%This code is desgined to simulate the Schlogl reaction system (6.1).
% This will reproduce figures 6.1 and 6.2 (a)
%By Kit Yates
%Created 19/10/15
%Last Modified 19/10/15

clear all
close all

%Define the rate constants
k1onu2=2.5*10^(-4); %(min^{-1})
k2onu=0.18; %(min^{-1})
k3=37.5; %(min^{-1})
k4nu=2200; %(min^{-1})

% Define the final time we will simulate to
T_final=100;

    % Define the number of repeats we will do
    M=10;
    
    
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
    
%Define the time after which we consider the process to be stationary
t_stat=2;

%Disregard the first t_stat seconds worth of data
rec_vector_cut=rec_vector_A(t_stat/rec_step:end,:);

%Find the maximum data value
max_data=max(rec_vector_cut(:));

%Decide how coarsely to plot the histograms
X_step=10;

%extract histogram data from each of the repeats after t_stat seconds when the
%process should have roughly reached steady state
[N,X]=hist(rec_vector_cut(:),[0:10:max_data]);

%Normalised the data to make it a PDF
N=N/(sum(N)*X_step);

% Plot the histogram of the stationary frequencies
[h1]=bar(X,N);

%correct the axes
axis([0 max_data+0.5 0 ceil(max(N)*100)/100])

%Define phi_0 to be one initially
phi(1)=1;


for i=1:max_data
    %Now find the frequencies for a range of values
    phi(i+1)=phi(i)*(k2onu*(i-1)*(i-2)+k4nu)/((k1onu2)*(i)*(i-1)*(i-2)+k3*i);
end

%Normalise the values of phi
phi=phi/sum(phi);

%Define the vector of determinsitc frequencies to plot
X_det=[0:1:max_data];

hold on
%Also plot the analytically determined mean
[h2]=plot(X_det,phi,'r','linewidth',5);
axis([0 max_data 0 0.015])

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('number of particles')
ylabel('stationary distribution')

legend([h1,h2],['Long-run SSA'],'Master equation');

exportfig(gcf,...
            ['schlogl_SPD.eps'],...
            'Format','eps2',...
            'Width','20',...
            'Color','cmyk',...
            'Resolution',300,...
            'FontMode','fixed',...
            'FontSize',21);
%             'LineWidth',6);
        %Save as a .fig as well
        saveas(gcf,['schlogl_SPD.fig'],'fig');