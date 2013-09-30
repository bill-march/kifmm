classdef DeltaKernel < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        tol;
        
    end
    
    methods
        
        function obj = DeltaKernel(eps)
            
            obj.tol = eps;
            
        end
        
        function res = eval(this, point1, point2)
           
            if (norm(point1-point2) < this.tol) 
               res = 1;
            else 
                res = 0;
            end
            
        end
        
        function res = eval_mat(this, data, rows, cols)
           
            res = zeros(size(rows,2), size(cols,2));
            
            for i = 1:size(rows,2)
               
                point_i = data(rows(i));
                
                for j = 1:size(cols,2)
                   
                    point_j = data(cols(j));
                    
                    res(i,j) = this.eval(point_i,point_j);
                    
                end
                
            end
            
        end
        
        
    end
    
end

