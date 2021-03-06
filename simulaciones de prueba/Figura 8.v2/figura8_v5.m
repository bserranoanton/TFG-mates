%This code is desgined to simulate the ...
%By Bel�n Serrano Ant�n
%Created 25/02/2020
%Last Modified --/--/2020

%Definici�n de par�metros:
t_cycle=0.15;     %Duraci�n del ciclo celular (desde el punto de
%restricci�n hasta la divisi�n
%celular)
t_apo=0.20;       %Duraci�n de la fase de ap�ptosis
t_next=0.3;       %Tiempo de interacci�n considerado en el modelo
alpha=10;         %Par�metro positivo dependiente del ant�geno (eq 5)
beta=0.1;         %Par�metro positivo dependiente del ant�geno (eq 5) 0.1

%Parameters for effector T cells
lambda_pd=0.04;    %Change in membrane receptor Rd, due to Rp signals 1 antes
lambda_taup=6*10^(-5); %Change in membrane receptor Rd, due to TCR signals
lambda_pp=0.5*10^(-5);%0.5*10^(-5);
mu_pc=0.4; %???????????? antes 0,4
mu_da=3.5; %??????????? antes 5

%Parameters for memory T cells
lambda_pd_mem = 0;  %Change in membrane receptor Rd, due to Rp signals
lambda_taup_mem = 10^(-6); %Change in membrane receptor Rd, due to TCR signals
lambda_pp_mem=2*10^(-2);%0.5*10^(-4);
mu_pc_mem=0.3; %???????????? 0.4

% Define the final time we will simulate to
T_final=25;

%Define the initial number of particles
N_init = 25; %tras la primera division habr� N_init effector y N_init memory
Y_init = 5;

%Define the initial number of T cell receptors
rp_init = 0;    %p0 no se usan porque ya tengo la matriz inicializada a 0
rd_init = 0;    %d0 

%Define how long the recording vectors should be
num_rec_steps=round(T_final/t_next);

%Initialise the index which will tell us where to write the current values
rec_ind=1;

%Define the maximum number of t cells
num_max_cells=10^7;

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
a0 = 0.3;
c0 = 0.08;% 0,01 nada mas salir
c0_mem = 0.04;
p0_mem = 0;
t_cell_matrix(1:2:2*N_init,1)=1;
t_cell_matrix(1:2:2*N_init,2)=0.3;%0.2; %a0
t_cell_matrix(1:2:2*N_init,3)=0.08;     %c0

t_cell_matrix(2:2:2*N_init,1)=2;
t_cell_matrix(2:2:2*N_init,3)=0.04; %c0 Hace que se disparen las memory
t_cell_matrix(2:2:2*N_init,4)=0;%p0 0.3;

% Initialise a vector which will hold the times when reactions occur
time_vec=zeros(1,num_rec_steps);

%Initialise the number of particles for this repeat
N_eff=N_init;
N_mem=N_init;
N=N_eff + N_mem;
Y=Y_init;

%Initialise index for t_cell_matrix
rec_ind_tcell_matrix=N+1;

%Define the initial time to be zero
t=0;


eliminado = 0;

muertas = 0;
creadas = 0;

while t<T_final
    
    %Increase the recording index
    rec_ind=rec_ind+1;
    
    if(eliminado==0)
        %Calculate Y
        Y = Y_init*exp(t*(alpha - N_eff*beta));
        Y = max(Y,0);
        if(Y < 10^(-6))
            Y=0;
            eliminado=1;
        end
    end
    
    %Fate decision for each T cell
    nCell=1;
    ind_N = 1;
    
    while nCell < rec_ind_tcell_matrix
        %Calculate rho
        %rho=rho_left*rand();
        %rho_left=rho_left-rho;
        %rho = 1/N;
        
        v_rand = rand(N,1)/N; %vector of N random numbers
       
        if(t_cell_matrix(nCell,1) == 1 || t_cell_matrix(nCell,1) == 2)
            rho = v_rand(ind_N);
            ind_N = ind_N + 1;
        end
       
        %rho = rho_left * (0.2*rand()); %para que una c�lula tenga como mucho el 20% de antigeno
        %rho_left=rho_left-rho;
        
        %Calculate r_tau
        r_tau=rho*Y;
        
      
        
        %Solve sys 9 (effector)
        if(t_cell_matrix(nCell,1) == 1 || t_cell_matrix(nCell,1) == 3)
            if(t_cell_matrix(nCell,6) > 0)
                %A�n est� en proceso de divisi�n
                t_cell_matrix(nCell,6) = max(t_cell_matrix(nCell,6)-t_next,0);
                if(t_cell_matrix(nCell,6)==0 && t_cell_matrix(nCell,1) == 3) %Se completa la division
                    N_eff=N_eff+1;
                    creadas = creadas +1;
                    t_cell_matrix(nCell,1) = 1;
                    disp('Hay division');
                end
            else
                %disp('Entramos a resolver sist 9');
                %C�lula lista para iniciar divisi�n (o ap�ptosis)
                %Condiciones iniciales
                p0_sys=t_cell_matrix(nCell,4);
                d0_sys=t_cell_matrix(nCell,5);
                c0_sys=t_cell_matrix(nCell,3);
                a0_sys=t_cell_matrix(nCell,2);
                
                %Calculamos c,a,p y d (Soluciones expl�citas de
                %sys9_sol(lambda_taup,lambda_pp, r_tau, p0, lambda_pd, d0, mu_pc, c0, mu_da, a0)
                [c,a,p,d] = sys9_sol(t,lambda_taup,lambda_pp, r_tau, p0_sys, lambda_pd, d0_sys, mu_pc, c0_sys, mu_da, a0_sys);
                
                if( a > 0 && c > 0)
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
                        
                        t_cell_matrix(nCell,3) = c0;
                        t_cell_matrix(nCell,2) = a0;
                        
                        t_cell_matrix(rec_ind_tcell_matrix,1)=3;
                        t_cell_matrix(rec_ind_tcell_matrix,4)=r_p_child_2;
                        t_cell_matrix(rec_ind_tcell_matrix,5)=r_d_child_2;
                        t_cell_matrix(rec_ind_tcell_matrix,6)=t_cycle;
                        
                        t_cell_matrix(rec_ind_tcell_matrix,2)=a0;
                        t_cell_matrix(rec_ind_tcell_matrix,3)=c0;
                        
                        rec_ind_tcell_matrix = rec_ind_tcell_matrix + 1;
                        
                        disp('Acabamos param div');
                    end
                end
            end
            nCell=nCell+1;
            
        elseif(t_cell_matrix(nCell,1) == 2 || t_cell_matrix(nCell,1) == 5) %memory
            if(t_cell_matrix(nCell,6) > 0)
                %A�n est� en proceso de divisi�n
                t_cell_matrix(nCell,6) = max(t_cell_matrix(nCell,6)-t_next,0);
                if(t_cell_matrix(nCell,6)==0 && t_cell_matrix(nCell,1) == 5) %Se completa la division
                    N_mem=N_mem+1;
                    t_cell_matrix(nCell,1) =2;
                    %disp('Hay division mem');
                end
            else
                %Soluciones de solutions_system_15.m
                p0_solsys =t_cell_matrix(nCell,4);
                c0_solsys =t_cell_matrix(nCell,3);
                
                [c,p] = sys15_sol(t,mu_pc_mem, p0_solsys, lambda_taup_mem, lambda_pp_mem, r_tau, c0_solsys);
                disp('P mem vale');
                disp(p);
                if(c <= 0)%dividimos
                    %disp('dividimos');
                    c=0;
                    delta_P_child_1 = 0.4+(0.6-0.4)*rand();
                    delta_P_child_2 = 1 - delta_P_child_1;
                   
                    r_p_child_1 = delta_P_child_1 * p;
                    r_p_child_2 = delta_P_child_2 * p;
                    
                    t_cell_matrix(nCell,4)=r_p_child_1;
                    t_cell_matrix(nCell,6)=t_cycle;
                    
                    t_cell_matrix(nCell,3) = c0_mem; %es 0
                    
                    t_cell_matrix(rec_ind_tcell_matrix,1)=5;
                    t_cell_matrix(rec_ind_tcell_matrix,4)=r_p_child_2;
                    t_cell_matrix(rec_ind_tcell_matrix,6)=t_cycle;
                    
                    t_cell_matrix(rec_ind_tcell_matrix,3)=c0_mem;
                    
                    rec_ind_tcell_matrix = rec_ind_tcell_matrix + 1;
                else
                    t_cell_matrix(nCell,4)=p;
                    t_cell_matrix(nCell,3)=c;
      
                end
            end
            
            nCell=nCell+1;
            %disp('Esta es memory');
            
        elseif(t_cell_matrix(nCell,1) == 4) %muerta effector
            %disp('Hay muerte 1');
            %disp(nCell);
            if(t_cell_matrix(nCell,6) > 0)
                %disp('Hay muerte');
                %disp(t_cell_matrix(nCell,6)-t_next);
                %A�n est� en proceso de muerte celular
                t_cell_matrix(nCell,6) = max(t_cell_matrix(nCell,6)-t_next,0);
                if(t_cell_matrix(nCell,6)==0) %Muere
                    disp('Hay muerte');
                    muertas = muertas +1;
                    disp(t_cell_matrix(nCell,6));
                    N_eff=N_eff-1;
                end
            end
            nCell=nCell+1;
        else
            break;
        end
        
    end %Acaba while de las c�lulas
    %Update the time
    t=t+t_next;
    
    %Record the time and the numbers of molecules
    time_vec(rec_ind)=t;
    N = N_eff + N_mem;
    rec_vector_N_eff(rec_ind)=N_eff;
    rec_vector_N_mem(rec_ind)=N_mem;
    rec_vector_Y(rec_ind)=Y;
    
end

disp('Creadas');
disp(creadas);

disp('muertas');
disp(muertas);

f1=figure;
%f2=figure;
%f3=figure;
%f4=figure;

%Draw results
figure(f1)
[hA1]=semilogy(time_vec,rec_vector_N_eff,'b');
%[hA1]=plot(time_vec,rec_vector_N_eff,'b');
hold on
[hA2]=semilogy(time_vec,rec_vector_Y,'r');
%[hA2]=plot(time_vec,rec_vector_Y,'r');
hold on
[hA3] = semilogy(time_vec,rec_vector_N_mem,'g');
%[hA3]=plot(time_vec,rec_vector_N_mem,'g');%semilogy(time_vec,rec_vector_N_mem,'g');
legend([hA1,hA3,hA2],'Effector T cells','Memory T cells','Pathogen');

% figure(f2)
% [hA1]=plot(time_vec,rec_vector_N_eff,'b');%semilogy(time_vec,rec_vector_N_eff,'b');
% legend(hA1,'Effector T cells');
% 
% figure(f3)
% [hA2]=plot(time_vec,rec_vector_Y,'r');%semilogy(time_vec,rec_vector_Y,'r');
% legend(hA2,'Pathogen');

%figure(f4)
%[hA3]=plot(time_vec,rec_vector_N_mem,'g');%semilogy(time_vec,rec_vector_N_mem,'g');
%legend(hA3,'Memory T cells');