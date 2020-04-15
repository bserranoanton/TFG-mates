%This code is desgined to produce data for the neural network
%By Belen Serrano Anton
%Created 15/04/2020
%Last Modified 15/04/2020

%Parameters
a = 3;
b = 0.1;
k = 0.4;
lambda = 0.5;

%Data file
fileID = fopen('data_neural_network_csv.csv','w');

while (a <= 3)
    b = 0.1;
    while(b <= 0.1)
        %Result of system 5.1
        [max_P, max_T, t_max_P, t_max_T, t_min_P, t_min_T] =  macro_func_neural_network(a, b, k, lambda);
        if(max_P ~= -1) %Intolerance
            disp(max_P);
            disp(max_T);
            disp(t_max_P);
            disp(t_max_T);
            disp(t_min_P);
            disp(t_min_T);
            fprintf(fileID,'%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n',max_P, max_T, t_max_P, t_max_T, t_min_P, t_min_T);
        end
        b = b + 0.1;
    end
    a = a + 0.1;
end

%Close file
fclose(fileID);