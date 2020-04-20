%This code is desgined to produce data for the neural network
%By Belen Serrano Anton
%Created 15/04/2020
%Last Modified 15/04/2020

%Parameters
%a = 0.85;
a = 0.1;
b = 0.1;
%k = 0.1;
k = 0.1;
lambda = 0.1;

%Data file
fileID = fopen('data_neural_network_csv_02.csv','w');

while (a <= 3)
    b = 0.1;
    while(b <= 3)
        k = 0.1;
        while(k <= 2)
            lambda = 0.1;
            while(lambda <= 2)
                %Result of system 5.1
                [max_P, max_T, t_max_P, t_max_T, t_min_P, t_min_T] =  macro_func_neural_network(a, b, k, lambda);
                if(max_P ~= -1) %Intolerance
                    fprintf(fileID,'%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n',max_P, max_T, t_max_P, t_max_T, t_min_P, t_min_T, a, b, k, lambda);
                end
                lambda = lambda + 0.3;
            end
            k = k + 0.3;
        end
        b = b + 0.2;
    end
    disp(a);
    a = a + 0.2;
end

%Close file
fclose(fileID);

