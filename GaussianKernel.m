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
        
        function [res, num_evals] = eval_mat(this, data, rows, cols, filter)
           
            res = zeros(size(rows,2), size(cols,2));
                
            num_evals = 0;
            
            if (nargin == 4)
              
                for i = 1:size(rows,2)
               
                    point_i = data(rows(i));

                    for j = 1:size(cols,2)

                        point_j = data(cols(j));

                        res(i,j) = this.eval(point_i,point_j);
                        num_evals = num_evals + 1;
                        
                    end

                end
                 
            else
                % then we have a filter
                [row_inds,col_inds] = find(filter);

                for i = 1:numel(row_inds)

                    row_point = data(rows(row_inds(i)));
                    col_point = data(cols(col_inds(i)));

                    res(row_inds(i), col_inds(i)) = this.eval(row_point, col_point);

                end
                num_evals = numel(row_inds);
                
            end
        end
        
    end
    
    
    
end