classdef GaussianKernel < handle
   
    properties
        
        Bandwidth 
    end
    
    
    methods
       
        function obj = GaussianKernel(band)
            
            obj.Bandwidth = band;
            
        end
        
        function res = eval(this, point1, point2)
           
            dist_sqr = (point1 - point2)^2;
            res = exp(-1 * dist_sqr * this.Bandwidth);
            
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