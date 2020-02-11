%This code is desgined to simulate the ...
%By Belén Serrano Antón
%Created 07/02/2020
%Last Modified --/--/2020


%Definición de parámetros:
t_cycle=0.15;     %Duración del ciclo celular (desde el punto de 
                    %restricción hasta la división
                    %celular)
t_apo=0.20;       %Duración de la fase de apóptosis
t_next=0.1;       %Tiempo de interacción considerado en el modelo
alpha=10;         %Parámetro positivo dependiente del antígeno (eq 5)
beta=0.015;       %Parámetro positivo dependiente del antígeno (eq 5)


%Parameters for effector T cells
lambda_pd=0.5;    %Change in membrane receptor Rd, due to Rp signals
lambda_taup=6*10^(-5); %Change in membrane receptor Rd, due to TCR signals
lambda_pp=0.5*10^(-5);
mu_pc=0.2; %????????????
mu_da=0.2; %??????????? 

%Parameters for memory T cells
lambda_pd_mem = 0;  %Change in membrane receptor Rd, due to Rp signals
lambda_taup_mem = 10^(-5); %Change in membrane receptor Rd, due to TCR signals
lambda_pp_mem=0.5*10^(-5); 
mu_pc_mem=0.2; %????????????

% Define the number of repeats we will do
%M=1;

% Define the final time we will simulate to
T_final=80;

%Define the initial number of particles
N_init = 100;
Y_init = 150;

%Define the initial number of T cell receptors
rp_init = 0;
rd_init = 0;

%Define how long the recording vectors should be
%num_rec_steps=10^6;
num_rec_steps=800;

%Initialise the index which will tell us where to write the current values
rec_ind=1;

%Define the maximum number of t cells
num_max_cells=10^6;

%Instantiate a vector which will hold the time varying values of N and Y
rec_vector_N_eff=-ones(1,num_rec_steps);
rec_vector_N_mem=-ones(1,num_rec_steps);
rec_vector_Y=-ones(1,num_rec_steps);

%Write the initial condision to these vectors
rec_vector_N_eff(rec_ind)=N_init/2; %Asymetric division
rec_vector_N_mem(rec_ind)=N_init/2;
rec_vector_Y(rec_ind)=Y_init;

% Initialise a vector which will hold the times when reactions occur
time_vec=-ones(1,num_rec_steps);

%Write the initial time to this vector
time_vec(rec_ind)=0;

%Define the initial time to be zero
t=0;

%Initialise the times to help recording
t_after=0;

%Initialise the number of particles for this repeat
N=N_init;
N_eff=N_init/2;
N_mem=N_init/2;
Y=Y_init;

%Initialise with N0 naive t cells
for k = 1:N
   rec_vector_tcell_2(k) = TCell(1,0,0);
end

%Initialise index
rec_ind_tcell_vector=N+1;

%Run through a for loop for each of the time-steps
%Enter the while loop
while t<20%T_final
    
    %Increase the recording index
    rec_ind=rec_ind+1;
    
    %Calculate Y
    ySol = sol_eq_5(alpha, beta, N,Y_init);
    Y=ySol(t);
    
    %Fate decision for each T cell
    nCell=1;
    rho_left=1;
    while nCell < 200%N
       
       %Calculate rho
       rho=rho_left*rand();
       rho_left=rho_left-rho;
       
       %Calculate r_tau
       r_tau=rho*Y;
       
       %Solve sys 9 (effector)
       if(rec_vector_tcell_2(nCell).My_type == 1)
                       
            p0= rec_vector_tcell_2(nCell).R_d;
            d0= rec_vector_tcell_2(nCell).R_p;
            c0=1;
            a0=1;
            %lambda_taup,lambda_pp, r_tau, p0, lambda_pd, d0, mu_pc, c0, mu_da, a0
            [cSol,aSol,pSol,dSol] = sols_sys_9_func(lambda_taup,lambda_pp, r_tau, p0, lambda_pd, d0, mu_pc, c0, mu_da, a0);
            
            if(aSol(t)<=0) %Apoptosis
                rec_vector_tcell_2(nCell).My_type=0;
                N_eff=N_eff-1;
            elseif(cSol(t)<=0) %Division
                children=rec_vector_tcell_2(nCell).divide();
                rec_vector_tcell_2(nCell)=children(1);
                rec_vector_tcell_2(rec_ind_tcell_vector+1)=children(2);
                rec_ind_tcell_vector=rec_ind_tcell_vector+1;
                N_eff=N_eff+1;
            end
            nCell=nCell+1; 
       %Solve sys 15 (memory)
       elseif(rec_vector_tcell_2(nCell).My_type == 2)
              nCell=nCell+1;
       %naive
       elseif(rec_vector_tcell_2(nCell).My_type == 3)
              nCell=nCell+1;
       %dead
       else
           
       end    
       %rec_ind_tcell_vector actualizar cuando se añada una célula
    end
    %Update the time
    t=t+t_next;
    
    %Record the time and the numbers of molecules
    time_vec(rec_ind)=t;
    N = N_eff + N_mem;
    rec_vector_N_eff(rec_ind)=N_eff;
    rec_vector_N_mem(rec_ind)=N_mem;    
    rec_vector_Y(rec_ind)=Y;
    
end

figure
    %plot the stochastic trajectory for A
    plot(time_vec,rec_vector_N_eff,'r');
    hold on
    plot(time_vec,rec_vector_Y,'b');
