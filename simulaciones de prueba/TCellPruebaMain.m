%Prueba TCell

% myCelula = TCell(1,0.1,0.2);
% 
% type = myCelula.My_type;
% 
% disp(type);
% 
% children = myCelula.divide();
% 
% disp(children(1).R_d);
% disp(children(1).R_p);
%format long

%[c,a,p,d] = sols_sys_9_func(6*10^(-5),0.5*10^(-5),2, 0, 0.5, 0, 0.2, 0, 0.2, 0);

%disp(c(0));
%sympref('FloatingPointOutput','default');

%disp(vpa(c(1)));

%Define the maximum number of t cells
num_max_cells=10^6;

%Instantiate a vector which will hold the T cells
%rec_vector_tcell = zeros(1, 10, TCell);
%objarray(100,1) = TCell(3,0,0);

for k = 1:105
   objarray(k) = TCell(3,0,0);
end

disp(objarray(105).My_type);
disp(objarray(105).R_p);