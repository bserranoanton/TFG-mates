%This code is designed to simulate the direction of locusts of two types.
%The locusts change direction in response to meeting two locusts going in
%the opposite direction with rate r and also change direction randomly with
%rate 1-r.
%By Kit Yates
% Created 17/10/13
%Last modified 04/04/14

close all
clear all

%Define the name we will use to save the data
save_name='locusts_minimal_SPD';

%Define the total number of individuals
N=20;

%Augment the save name
save_name=[save_name,'N=',num2str(N)];

%Define the number of right moving individuals initially
X1=ceil(N*rand);

%Define the number of left moving individuals initially
X2=N-X1;

%***These rates are the avergaed rates fitted to the experimental data***

%Define the rate at which particles switch from one state to the other
%spontaneously
r1=0.0225;

%Augment the save name
save_name=[save_name,'r1=',num2str(r1)];

%Define the rates at which particles switch as a result of encounters
r2=0.0453;

%Augment the save name
save_name=[save_name,'r2=',num2str(r2)];

r3=0.1664;

%Augment the save name
save_name=[save_name,'r3=',num2str(r3)];

r4=0;

%Augment the save name
save_name=[save_name,'r4=',num2str(r4)];

%Define a final time for the simulation
T_final=100000000;

%Define the initial time to be zero
t=0;

%Define an appropriate time-step to record X1 at
rec_time_step=1;

%Define how many time_points it is possible to record at this frequency.
num_rec_steps=T_final/rec_time_step;

%initialise a vector to record the value of X1
%Define the vector which records the number of times each aligment has been
%visited
X1_rec=zeros(N+1,1);

%Enter a while loop that we don't leave until we get past the last time
while t<T_final
    
    %Define the propensity for a right moving individual to become a left
    %moving individual
    alr=(r1*X1/N+r2*X1*X2/N^2+(r3*X1*X2^2)/N^3+r4*X1*X2^3/N^4);
    arl=(r1*X2/N+r2*X1*X2/N^2+(r3*X2*X1^2)/N^3+r4*X2*X1^3/N^4);
    
    %Find the sum of the propensity functions
    suma=arl+alr;
    
    %Determine which reaction will take place
    if (rand*suma)<alr
        %Then convert a right moving individual to a left moving individual
        X1=X1-1;
        X2=X2+1;
    else
        %Then convert a left moving indicidual to a right-moving individual
        X2=X2-1;
        X1=X1+1;
    end
    
    %Calculate the time-step
    tau=log(1/rand)/suma;
    
    %Update the timestep
    t=t+tau;
    
    %Calculate the index of the recording vector at the start of the time-step
    %we add eps to ensure we never get a value of zero for min_ind
    min_ind=ceil((t-tau+eps)/rec_time_step);
    
    %Calculate the index of the recording vector at the end of the time-step
    max_ind=floor(t/rec_time_step);
    
    %record the cell density every rec_time_step units
    if min_ind<=max_ind
        
        %find the number of entries we are going to make in each of the
        %vectors
        num_entries=max_ind-min_ind+1;
        
        X1_rec(X1+1)=X1_rec(X1+1)+num_entries;
    end
    
end

%Add .mat to the save_name
save_name=[save_name,'.mat'];


%Save the data with the relevant name
save(save_name)