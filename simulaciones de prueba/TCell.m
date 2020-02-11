classdef TCell
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %Effector = 1
        %Memory = 2
        %Muerta = 0
        My_type
        %Delta_P
        %Delta_D
        %R_tau
        R_p
        R_d
    end
    
    methods
        function obj = TCell(type, r_p, r_d)
            obj.My_type = type;
            obj.R_p = r_p;
            obj.R_d = r_d;
        end
        
        function children = divide(mother)
            delta_P_child_1 = randi([4 6]) / 10;
            delta_P_child_2 = 1 - delta_P_child_1;
            delta_D_child_1 = randi([4 6]) / 10;
            delta_D_child_2 = 1 - delta_D_child_1;
            
            r_p_child_1 = delta_P_child_1 * mother.R_p;
            r_p_child_2 = delta_P_child_2 * mother.R_p;
            
            r_d_child_1 = delta_D_child_1 * mother.R_d;
            r_d_child_2 = delta_D_child_2 * mother.R_d;
            
            children(1) = TCell(mother.My_type,r_p_child_1,r_d_child_1);
            children(2) = TCell(mother.My_type,r_p_child_2,r_d_child_2);
        end
    end
    
end

