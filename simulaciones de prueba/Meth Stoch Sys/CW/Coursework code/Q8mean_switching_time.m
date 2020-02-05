%This code is desgined to calculate the mean switching time for the system
%question 8, coursework 2017.
%Candidate number: 19819

clear all
close all

%Define the rate constant of degradation
k1=2.5; %(sec^{-1})
k2onu=2; %(sec^{-1})
k3onu2=0.005; %(sec^{-1})

% Define the number of repeats we will do
M=100; 

%Define the number of repeats to record
M_rec=100;

%Define the number of particles
N = 50;

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
A_init=2;
B_init=48;


%Define the initial number of molecules forthe recording vector
rec_vector_A(1,:)=A_init*ones(1,M);
rec_vector_B(1,:)=B_init*ones(1,M);

%Define a vector which will hold the propensity functions
a=zeros(6,1);

%square periods so we can calculate the sample variance easily
sum_time_switching=0;
sum_time_switching_square=0;

%Initialise a counter to say we are on the first period
time_switching_counter=0;

%Initialise sample vairance to be very large so we enter the repeat while
%loop
var_sample_mean=inf;

%Initialise a repeat counting variable to be zero
i=0;

time_switching_length=zeros(1,M*T_final/2);

%Run through a for loop for each of the repeats
while var_sample_mean>10^(-6) || i<M_rec;
    
     %Increase the counting index
    i=i+1;
    %Define the initial time to be zero
    t=0;
    
    %Initialise the times to help recording
    t_after=0;
    
    %   initialise the number of particles for this repeat
    A=A_init;
    B=B_init;
    
    %Run through a for loop for each of the time-steps
    %Enter the while loop
    while t<T_final && A < N/2
        
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
    time_switching_counter = time_switching_counter+1;
    sum_time_switching = sum_time_switching + t;
    sum_time_switching_square = sum_time_switching_square + t^2;
    %Calculate the sample variance of the period
    var_sample_mean=(sum_time_switching_square/time_switching_counter-(sum_time_switching/time_switching_counter)^2)/time_switching_counter;
    
end
%Get rid of the first couple of entries of the switching time in each
%repeat
mean_time_switching=sum_time_switching/time_switching_counter;
disp('Mean time switching');
disp(mean_time_switching);