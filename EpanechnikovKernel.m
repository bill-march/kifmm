classdef EpanechnikovKernel < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        bandwidth
        
    end
    
    methods
        
        function obj = EpanechnikovKernel(band)
        
            obj.bandwidth = band;
            
        end
        
        function res = eval(this, point1, point2)
           
            distsqr = (point1 - point2)' * (point1 - point2) / (this.bandwidth * this.bandwidth);
            if (distsqr < 1)
               res = 3/4 * (1 - distsqr);
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

