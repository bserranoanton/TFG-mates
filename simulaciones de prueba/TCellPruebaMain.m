%Prueba TCell

myCelula = TCell(1,0.1,0.2);

type = myCelula.My_type;

disp(type);

children = myCelula.divide();

disp(children(1).R_d);
disp(children(1).R_p);