%define the rate  constant of degradation
k = 0.1; %sec^-1

%def number of repeats
M=100;
%define the final time we will simulate to
T_final =30;

%define the recording time interval
rec_step =0.01;

%define avector of time points
time_vec = [0:rec_step:T_final];

%Calculate the number of recording steps required
num_rec_steps = T_final/rec_step;

%matrix witch will record the number of particles each time point at ecah
%time at each repeat
rec_vector = zeros(num_rec_steps+1,M);

%initial number of particles
A_init=20;

%initiall numer of molecules for recording vector
rec_vector(1,:) = A_init*ones(1,M);


%run throug a for loop for each of the repeats

for i= 1:M
    %define the initial time to be zero
    t=0;
    
    %initialise the times to help recording
    t_after =0;
    
    %initialise the number of particles for this repeat
    A = A_init;
    
    %enter while loop
    while t < T_final
        %calculate the propensity functions
        a0=A*k;
        
        %determinate the time for the next reaction
        tau = (1/a0)*log(1/rand);
        
        %update the time
        t = t+tau;
        
        %write the current value of A to A_old
        A_old=A;
        
        %implement the degradation reaction
        A = max(A_old-1,0);
        
        %calculate the times for recording
        t_before = t_after;
        t_after =t;
        
        %calculate the indices of the time step before and the current time
        %step in terms of recording
        ind_before = ceil((t_before+eps)/rec_step); %eps => nozero
        ind_after = min(floor(t_after/rec_step),num_rec_steps);
        
        %find out how many time-steps to write to
        steps_to_write = ind_after - ind_before+1; 
        
        if steps_to_write > 0 && steps_to_write ~=Inf
        rec_vector(ind_before+1:ind_after+1,i) = A_old*ones(steps_to_write,1);
        end
        
        %End the while loop
    end
    
    %End the for loop
end


%Calculate the mean over all 100 realisations
mean_A = mean(rec_vector,2);

%calculate the analytically determines mean from the master eq
analytical_mean = A_init*exp(-k*time_vec);

%Plot the first 10 trajectories


%Plot the first 10 trajectories
if M>=10
    plot(time_vec,rec_vector(:,1:10));
else
    plot(time_vec,rec_vector);
end

hold on
%Also plot on the mean over 100 realisations
[h1]=plot(time_vec,mean_A,'r','linewidth',5);

%Also plot the analytically determined mean
[h2]=plot(time_vec,analytical_mean,'k--','linewidth',5);

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of particles')

legend([h1,h2],['Mean of M=',num2str(M),' simulations'],'Analytical mean');

exportfig(gcf,...
    ['degradation_realisations_and_mean.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',21);
%             'LineWidth',6);
%Save as a .fig as well
saveas(gcf,['degradation_realisations_and_mean.fig'],'fig');

%Plot the figures from lecture 3 using the same code
%Make a new figure
figure

%Also plot the analytically determined mean
plot(time_vec,analytical_mean,'k--','linewidth',5);

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of particles')
exportfig(gcf,...
    ['production_degradation_mean.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',21);

%Make a new figure
figure

%plot the analytically determined mean
[h2]=plot(time_vec,analytical_mean,'k--','linewidth',5);

legend(h2,'Analytical mean');


%Plot the first trajectory
plot(time_vec,rec_vector(:,1));

% %Set the title
% title(['t= ', num2str(t*rec_step_increment)])
%Set the x and y labels
xlabel('time (sec)')
ylabel('number of particles')
exportfig(gcf,...
    ['production_degradation_mean_and_trajectory.eps'],...
    'Format','eps2',...
    'Width','20',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',21);






