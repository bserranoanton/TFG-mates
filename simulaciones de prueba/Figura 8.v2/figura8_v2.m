%This code is desgined to simulate the ...
%By Belén Serrano Antón
%Created 11/02/2020
%Last Modified --/--/2020


%Definición de parámetros:
t_cycle=0.15;     %Duración del ciclo celular (desde el punto de
%restricción hasta la división
%celular)
t_apo=0.20;       %Duración de la fase de apóptosis
t_next=0.3;       %Tiempo de interacción considerado en el modelo
alpha=10;         %Parámetro positivo dependiente del antígeno (eq 5)
beta=0.3;         %Parámetro positivo dependiente del antígeno (eq 5)

%Parameters for effector T cells
lambda_pd=0.5;    %Change in membrane receptor Rd, due to Rp signals
lambda_taup=6*10^(-5); %Change in membrane receptor Rd, due to TCR signals
lambda_pp=0.5*10^(-5);
mu_pc=0.2; %????????????
mu_da=3; %???????????

%Parameters for memory T cells
lambda_pd_mem = 0;  %Change in membrane receptor Rd, due to Rp signals
lambda_taup_mem = 10^(-5); %Change in membrane receptor Rd, due to TCR signals
lambda_pp_mem=0.5*10^(-5);
mu_pc_mem=2; %????????????

% Define the number of repeats we will do
M=1;

% Define the final time we will simulate to
T_final=100;

%Define the initial number of particles
N_init = 25; %tras la primera division habrá N_init effector y N_init memory
Y_init = 5;

%Define the initial number of T cell receptors
rp_init = 0;
rd_init = 0;

%Define how long the recording vectors should be
num_rec_steps=round(T_final/t_next);

%Initialise the index which will tell us where to write the current values
rec_ind=1;

%Define the maximum number of t cells
num_max_cells=10^6;

%Instantiate a vector which will hold the time varying values of N and Y
rec_vector_N_eff=-ones(1,num_rec_steps);
rec_vector_N_mem=-ones(1,num_rec_steps);
rec_vector_Y=-ones(1,num_rec_steps);

%Write the initial condision to these vectors
rec_vector_N_eff(rec_ind)=N_init; %Asymetric division
rec_vector_N_mem(rec_ind)=N_init;
rec_vector_Y(rec_ind)=Y_init;

%Instantiate a vector which will hold the t cells
t_cell_matrix=zeros(num_max_cells,6);

%Write the initial condision to this vector
t_cell_matrix(1:2:2*N_init,1)=1;
t_cell_matrix(1:2:2*N_init,2)=10;
t_cell_matrix(1:2:2*N_init,3)=1;

t_cell_matrix(2:2:2*N_init,1)=2;
t_cell_matrix(2:2:2*N_init,3)=1;

%disp(t_cell_matrix(1:N_init,1:4));

% Initialise a vector which will hold the times when reactions occur
time_vec=-ones(1,num_rec_steps);


f1=figure;
%f2=figure;
for repeats=1:M
    
    %Hay código repetido de arriba para resetear valores, habría que
    %ponerlo solo aquí.
    
    %Write the initial time to this vector
    time_vec(rec_ind)=0;
    
    %Define the initial time to be zero
    t=0;
    
    %Initialise the times to help recording
    t_after=0;
    
    %Initialise the number of particles for this repeat
    N_eff=N_init;
    N_mem=N_init;
    N=N_eff + N_mem;
    Y=Y_init;
    
    %Initialise index for t_cell_matrix
    rec_ind_tcell_matrix=N+1;
    
    
    %Write the initial condision to these vectors
    rec_vector_N_eff(rec_ind)=N_eff; %Asymetric division
    rec_vector_N_mem(rec_ind)=N_mem;
    rec_vector_Y(rec_ind)=Y_init;
    
    %Write the initial condision to this vector
    t_cell_matrix(1:2:2*N_init,1)=1;
    t_cell_matrix(1:2:2*N_init,2)=10;
    t_cell_matrix(1:2:2*N_init,3)=1;

    t_cell_matrix(2:2:2*N_init,1)=2;
    t_cell_matrix(1:2:2*N_init,3)=1;
    
    rec_ind=1;
    
    
    %Run through a for loop for each of the time-steps
    %Enter the while loop
    while t<T_final
        
        %Increase the recording index
        rec_ind=rec_ind+1;
        
        %Calculate Y
        Y = Y_init*exp(t*(alpha - N_eff*beta));
        Y = max(Y,0);
        
        %Fate decision for each T cell
        nCell=1;
        rho_left=1;
        %disp('**********************************');
        while nCell < N
            
            %disp('Vamos por la celula ');
            %disp(nCell);
            %Calculate rho
            rho=rho_left*rand();
            rho_left=rho_left-rho;
            
            %rho = 1/N;
            %disp('Acaba de calcular rho');
            %Calculate r_tau
            r_tau=rho*Y;
            
            %Solve sys 9 (effector)
            if(t_cell_matrix(nCell,1) == 1)
                %disp('Efectora');
                if(t_cell_matrix(nCell,6) > 0)
                    %Aún está en proceso de división
                    t_cell_matrix(nCell,6) = max(t_cell_matrix(nCell,6)-t_next,0);
                    if(t_cell_matrix(nCell,6)==0) %Se completa la division
                        N_eff=N_eff+1;
                        disp('Hay division');
                    end
                else
                    %disp('Entramos a resolver sist 9');
                    %Célula lista para iniciar división (o apóptosis)
                    %Condiciones iniciales
                    p0=t_cell_matrix(nCell,4);
                    d0=t_cell_matrix(nCell,5);
                    c0=t_cell_matrix(nCell,3);
                    a0=t_cell_matrix(nCell,2);
                    
                    %Calculamos c,a,p y d (Soluciones explícitas de
                    %solutions_system9.m)
                    p=(lambda_taup*r_tau + exp(-lambda_pp*t)*(lambda_pp*p0 - lambda_taup*r_tau))/lambda_pp;
                    d=d0 + (lambda_pd*p0 - (lambda_pd*lambda_taup*r_tau)/lambda_pp)/lambda_pp - ...
                        (exp(-lambda_pp*t)*(lambda_pd*p0 - (lambda_pd*lambda_taup*r_tau)/lambda_pp) - ...
                        lambda_pd*lambda_taup*r_tau*t)/lambda_pp;
                    c=c0 - (mu_pc*p0 - (lambda_taup*mu_pc*r_tau)/lambda_pp)/lambda_pp + ...
                        (exp(-lambda_pp*t)*(mu_pc*p0 - (lambda_taup*mu_pc*r_tau)/lambda_pp) - ...
                        lambda_taup*mu_pc*r_tau*t)/lambda_pp;
                    a=a0 + (lambda_pd*mu_da*p0 - (lambda_pd*lambda_taup*mu_da*r_tau)/lambda_pp)/lambda_pp^2 ...
                        - (t*(d0*mu_da*lambda_pp^2 + lambda_pd*mu_da*p0*lambda_pp - ...
                        lambda_pd*lambda_taup*mu_da*r_tau) + exp(-lambda_pp*t)*(lambda_pd*mu_da*p0 ...
                        - (lambda_pd*lambda_taup*mu_da*r_tau)/lambda_pp) + ...
                        (lambda_pd*lambda_pp*lambda_taup*mu_da*r_tau*t^2)/2)/lambda_pp^2;
                    
                    if((d <=0 || p<=0) && a >= 0 && c >= 0)
                        d=max(d,0);
                        p=max(p,0);
                        t_cell_matrix(nCell,4) = p;
                        t_cell_matrix(nCell,5) = d;
                        t_cell_matrix(nCell,3) = c;
                        t_cell_matrix(nCell,2) = a;
                    else
                        if(a <= 0) %apoptosis
                            a=0;
                            disp('Matamos');
                            disp(nCell);
                            t_cell_matrix(nCell,6)=t_apo;
                            t_cell_matrix(nCell,1)=4;
                        elseif(c <= 0 && a > 0) %division
                            %si dejo condición a > 0 se dispara el patógeno
                            disp('dividimos Effector');
                            c=0;
                            %delta_P_child_1 = randi([4 6]) / 10;
                            %Random number between 0.4 and 0.6
                            delta_P_child_1 = 0.4+(0.6-0.4)*rand();
                            delta_P_child_2 = 1 - delta_P_child_1;
                            delta_D_child_1 = 0.4+(0.6-0.4)*rand();
                            delta_D_child_2 = 1 - delta_D_child_1;
                            
                            r_p_child_1 = delta_P_child_1 * p;
                            r_p_child_2 = delta_P_child_2 * p;
                            
                            r_d_child_1 = delta_D_child_1 * d;
                            r_d_child_2 = delta_D_child_2 * d;
                            
                            t_cell_matrix(nCell,4)=r_p_child_1;
                            t_cell_matrix(nCell,5)=r_d_child_1;
                            t_cell_matrix(nCell,6)=t_cycle;
                            
                            t_cell_matrix(nCell,3) = c;
                            t_cell_matrix(nCell,2) = a;
                            
                            t_cell_matrix(rec_ind_tcell_matrix,1)=1;
                            t_cell_matrix(rec_ind_tcell_matrix,4)=r_p_child_2;
                            t_cell_matrix(rec_ind_tcell_matrix,5)=r_d_child_2;
                            t_cell_matrix(rec_ind_tcell_matrix,6)=t_cycle;
                            
                            t_cell_matrix(rec_ind_tcell_matrix,2)=a;
                            t_cell_matrix(rec_ind_tcell_matrix,3)=c;
                            
                            rec_ind_tcell_matrix=rec_ind_tcell_matrix+1;
                            
                            disp('Acabamos param div');
                        end
                    end
                end
                nCell=nCell+1;
                
            elseif(t_cell_matrix(nCell,1) == 2) %memory
                if(t_cell_matrix(nCell,6) > 0)
                    %Aún está en proceso de división
                    t_cell_matrix(nCell,6) = max(t_cell_matrix(nCell,6)-t_next,0);
                    if(t_cell_matrix(nCell,6)==0) %Se completa la division
                        N_mem=N_mem+1;
                        %disp('Hay division mem');
                    end
                else
                    %Soluciones de solutions_system_15.m
                    p0=t_cell_matrix(nCell,4);
                    c0=t_cell_matrix(nCell,3);
                    c=c0 - (mu_pc*p0 - (lambda_taup*mu_pc*r_tau)/lambda_pp)/lambda_pp + (exp(-lambda_pp*t)*(mu_pc*p0 - (lambda_taup*mu_pc*r_tau)/lambda_pp) - lambda_taup*mu_pc*r_tau*t)/lambda_pp;
                    p=(lambda_taup*r_tau + exp(-lambda_pp*t)*(lambda_pp*p0 - lambda_taup*r_tau))/lambda_pp;
                    
                    if(c <= 0)%dividimos
                        %disp('dividimos');
                        c=0;
                        delta_P_child_1 = 0.4+(0.6-0.4)*rand();
                        delta_P_child_2 = 1 - delta_P_child_1;
                        %delta_D_child_1 = 0.4+(0.6-0.4)*rand();
                        %delta_D_child_2 = 1 - delta_D_child_1;
                        
                        r_p_child_1 = delta_P_child_1 * p;
                        r_p_child_2 = delta_P_child_2 * p;
                        
                        %r_d_child_1 = delta_D_child_1 * d;
                        %r_d_child_2 = delta_D_child_2 * d;
                        
                        t_cell_matrix(nCell,4)=r_p_child_1;
                        %t_cell_matrix(nCell,5)=r_d_child_1;
                        t_cell_matrix(nCell,6)=t_cycle;
                        
                        t_cell_matrix(nCell,3) = c;
                        
                        t_cell_matrix(rec_ind_tcell_matrix,1)=2;
                        t_cell_matrix(rec_ind_tcell_matrix,4)=r_p_child_2;
                        %t_cell_matrix(rec_ind_tcell_matrix,5)=r_d_child_2;
                        t_cell_matrix(rec_ind_tcell_matrix,6)=t_cycle;
                        
                        t_cell_matrix(rec_ind_tcell_matrix,3)=c;
                        
                        rec_ind_tcell_matrix=rec_ind_tcell_matrix+1;
                    end
                end
                
                nCell=nCell+1;
                %disp('Esta es memory');
                
            elseif(t_cell_matrix(nCell,1) == 4) %muerta effector
                %disp('Hay muerte 1');
                %disp(nCell);
                if(t_cell_matrix(nCell,6) > 0)
                    disp('Hay muerte');
                    disp(t_cell_matrix(nCell,6)-t_next);
                    %Aún está en proceso de muerte celular
                    t_cell_matrix(nCell,6) = max(t_cell_matrix(nCell,6)-t_next,0);
                    if(t_cell_matrix(nCell,6)==0) %Muere
                        N_eff=N_eff-1;
                        
                    end
                end
                nCell=nCell+1;
            else
                break;
            end
            
        end %Acaba while de las células
        %Update the time
        t=t+t_next;
        
        %Record the time and the numbers of molecules
        time_vec(rec_ind)=t;
        N = N_eff + N_mem;
        rec_vector_N_eff(rec_ind)=N_eff;
        rec_vector_N_mem(rec_ind)=N_mem;
        rec_vector_Y(rec_ind)=Y;
        
    end
    
    %Draw results
    figure(f1)
    [hA1]=semilogy(time_vec,rec_vector_N_eff,'b');
    hold on
    %figure(f2)
    [hA2]=semilogy(time_vec,rec_vector_Y,'r');
    %plot(time_vec,rec_vector_Y,'b');
    %hold on
    %[hA3]=semilogy(time_vec,rec_vector_N_mem,'g');
end

%figure(f1)
%legend([hA1,hA3,hA2],'Effector T cells','Memory T cells','Pathogen');
%plot the stochastic trajectory for A
%plot(time_vec,rec_vector_N_eff,'r');

