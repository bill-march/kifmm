classdef UnitaryKernel < handle
    % kernel where all interactions are 1
    
    properties
    end
    
    methods

        function obj = UnitaryKernel()
            
        end
        
        function res = eval(this, point1, point2)
           
            res = 1.0;
            
        end
        
        function res = eval_mat(this, data, rows, cols, filter)
           
            res = ones(size(rows,2), size(cols,2));
            if (nargin > 4)
                % this assumes the values in filter are ones
                res = res .* filter;
            end
        end
    
    end
    
end

