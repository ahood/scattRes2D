classdef schurComp < handle
    properties
        I1
        I2
    end
    methods
        function obj = schurComp(I1,I2)
            obj.I1 = I1; obj.I2 = I2; 
        end
        function s = S(obj,X)
            I1 = obj.I1; I2 = obj.I2;   
            s = X(I1,I1)-X(I1,I2)*(X(I2,I2)\X(I2,I1));
        end 
        function u = update(obj,X)
            I1 = obj.I1; I2 = obj.I2;
            u = -X(I1,I2)*(X(I2,I2)\X(I2,I1));
        end
        function s = SRHS(obj,X,RHS)
            I1 = obj.I1; I2 = obj.I2;
            s = RHS(I1) - X(I1,I2)*(X(I2,I2)\RHS(I2));
        end
    end
end