%This code is desgined to plot the stationary probability distribution of the system,
%as found by long run stochastic simulation, as a histogram,
%and the analytically derived stationary probability distribution. 
%Question 9, coursework 2017
%Candidate number: 19819

clear all
close all

%Define the rate constant of degradation
k1=2.5; %(sec^{-1})
k2onu=2; %(sec^{-1})
k3onu2=0.003; %(sec^{-1})

% Define the number of repeats we will do
M=1000; 

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
A_init=25;
B_init=25;


%Define the initial number of molecules forthe recording vector
rec_vector_A(1,:)=A_init*ones(1,M);
rec_vector_B(1,:)=B_init*ones(1,M);

%Define a vector which will hold the propensity functions
a=zeros(6,1);

%Run through a for loop for each of the repeats
for i=1:M
    
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
    
end

%Define the time after which we consider the process to be stationary
t_stat=20;

%Disregard the first t_stat seconds worth of data
rec_vector_A_cut=rec_vector_A(t_stat/rec_step:end,:);
rec_vector_B_cut=rec_vector_B(t_stat/rec_step:end,:);

%Calculate the maximum values of the A and B data
max_data_A=max(rec_vector_A_cut(:));
max_data_B=max(rec_vector_B_cut(:));

%Decide the range into which the data will fall
X_A=[0:1:max_data_A];
X_B=[0:1:max_data_B];

%Create a matrix which will hold the frequencies of occurance in each
%pairing of A and B
N=zeros(max_data_A+1,max_data_B+1);

%Create the histogram data manually
for i=1:length(rec_vector_A_cut)
    for j=1:M
        N(rec_vector_A_cut(i,j)+1,rec_vector_B_cut(i,j)+1)=N(rec_vector_A_cut(i,j)+1,rec_vector_B_cut(i,j)+1)+1;
    end
end
%Normalised the data to make it a PDF
N=N/sum(N(:));

%Decide how much of the data to show
max_show_data=50;

%Define the values of the data to plot
X=[0:1:max_show_data];

%extract the parts of the histogram to show on the plot
N_show=N(1:max_show_data+1,1:max_show_data+1);

% Plot the histogram of the stationary frequencies
[h1]=pcolor(X,X,N_show');
%Note pcolor isn't exactly plotting the histogram here, but interpolating
%between points, but it is close enough to give the feel for it.

Z = 50;
%calculate SPD analytically
 phi = zeros(1,51);
    phi(1) = 1;
    phi(2) = phi(1)*((k1*Z)/(k1+k2onu*(Z-1)+k3onu2*(Z-1)*(Z-2)));
   
    %recursion
    for n = 2:50
        coeffn = k1*Z + 2*k2onu*n*(Z-n) + k3onu2*n*(Z-n)*(Z-2);
        coeffnminus1 = k1*(Z-n+1)+k2onu*(n-1)*(Z-n+1)+k3onu2*(n-1)*(n-2)*(Z-n+1); 
        coeffnplus1 = k1*(n+1)+k2onu*(n+1)*(Z-n-1)+k3onu2*(n+1)*(Z-n-1)*(Z-n-2);
        phi(n+1)=  (phi(n)*coeffn - phi(n-1)*coeffnminus1)/coeffnplus1;
    end
    %normalise
   phi=phi/sum(phi);


%colorbar

%Change the shading
% shading interp

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('number of A particles')
ylabel('number of B particles')

% exportfig(gcf,...
%             ['dimerisation_SPD.eps'],...
%             'Format','eps2',...
%             'Width','20',...
%             'Color','cmyk',...
%             'Resolution',300,...
%             'FontMode','fixed',...
%             'FontSize',21);
% %             'LineWidth',6);
%         %Save as a .fig as well
%         saveas(gcf,['dimerisation_SPD.fig'],'fig');

%Now calculate the marginal stationary probability distribution for A
N_A=sum(N,2);

%Truncate N_A to just show the relevant range of values
N_A_show=N_A(1:max_show_data+1);

%Make a new figure
figure

%plot the histogram for the SPD
[h2]=bar(X,N_A_show);

%correct the axes
axis([0 max_show_data+0.5 0 ceil(max(N_A_show)*100)/100])

%Decide which value of the histogram we want to highlight
%patch_highlight=10;

hold on
%plot a patch over the bar which corresponds to the mean-field mean value
%for A
%patch([patch_highlight-0.41 patch_highlight+0.41 patch_highlight+0.41 patch_highlight-0.41],[0 0 N_A_show(patch_highlight+1) N_A_show(patch_highlight+1)],'r')
plot(phi, 'linewidth',5);
% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('number of A particles')
ylabel('stationary distribution')

% exportfig(gcf,...
%             ['dimerisation_marginal_SPD_A.eps'],...
%             'Format','eps2',...
%             'Width','20',...
%             'Color','cmyk',...
%             'Resolution',300,...
%             'FontMode','fixed',...
%             'FontSize',21);
%             'LineWidth',6);
%         %Save as a .fig as well
%         saveas(gcf,['dimerisation_marginal_SPD_A.fig'],'fig');