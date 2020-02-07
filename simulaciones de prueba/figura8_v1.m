%This code is desgined to simulate the ...
%By Belén Serrano Antón
%Created 07/02/2020
%Last Modified --/--/2020


%Define params:
t_cycle = 0.15;
t_apo = 0.20;
alpha = 10;
beta = 0.015;

%For effector T cells
lambda_pd = 0.5;
gamma_lambda_taup = 6*10^(-5);

%For memory T cells
lambda_pd_mem = 0;
gamma_lambda_taup_mem = 10^(-5);

% Define the number of repeats we will do
%M=1;

% Define the final time we will simulate to
T_final=80;

%Define the initial number of particles
N_init = 100;
Y_init = 200;

%Define the initial number of T cell receptors
rp_init = 0;
rd_init = 0;

%Define how long the recording vectors should be
num_rec_steps=10^6;

%Initialise the index which will tell us where to write the current values
rec_ind=1;

%Instantiate a vector which will hold the time varying values of N and Y
rec_vector_N=-ones(1,num_rec_steps);
rec_vector_Y=-ones(1,num_rec_steps);

%Write the initial condisiotn to these vectors
rec_vector_N(rec_ind)=N_init;
rec_vector_Y(rec_ind)=Y_init;

%Instantiate a vector which will hold the time varying values the number of
%Ri receptors in T cell membrane
rec_rp_vector_N=-ones(1,num_rec_steps);
rec_rd_vector_N=-ones(1,num_rec_steps);

%Write the initial condisiotn to these vectors
rec_rp_vector_N(rec_ind)=rp_init;
rec_rd_vector_N(rec_ind)=rd_init;

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
N=A_init;
Y=B_init;

%Run through a for loop for each of the time-steps
%Enter the while loop
while t<T_final
    
    %Increase the recording index
    rec_ind=rec_ind+1;
    
    %Define the propensity functions
    a(1)=A*(A-1)*B*k1onu2;
    a(2)=k2nu;
    a(3)=k3*A;
    a(4)=k4nu;
    
    %Calculate the cumulative sum of a
    cumsuma=cumsum(a);
    
    %Calculate the propensity functions
    a0=cumsuma(end);
    
    %Determine the time for the next reaction
    t_next= %Depende de si se divide o muere
    
    %Update the time
    t=t+t_next;
    
    %Write the current value of A to A_old
    N_old=N;
    Y_old=Y;
    
    %multiply the random number by the sum of the propensities and
    %store so we can reuse
    ra0=rand*a0;
    
    %if we have not run out of molecules
    if t_next<inf
        %Determine which reaction happened
        if ra0<cumsuma(1)
            %Implement the thrd order production reaction
            A=A_old+1;
            B=B_old-1;
        elseif ra0<cumsuma(2)
            %Implement the zeroth order production reaction
            A=A_old+1;
        elseif ra0<cumsuma(3)
            % Implement a first order degradation reaction
            A=A_old-1;
        else
            % Implement the zeroth order poduction reaction for B
            B=B_old+1;
        end
    end
    
    %Rcord the time and the numbers of molecules
    time_vec(rec_ind)=t;
    rec_vector_N(rec_ind)=N;
    rec_vector_Y(rec_ind)=Y;
    
end

