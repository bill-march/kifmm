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
        
    end
    
    
    
end