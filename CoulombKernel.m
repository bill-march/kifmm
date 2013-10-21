classdef CoulombKernel < handle
   
    properties
        
        Scale 
        
    end
    
    
    methods
       
        function obj = CoulombKernel(scale)
            
            if (nargin > 0)
                obj.Scale = scale;
            else 
                obj.Scale = 1;
            end
            
        end
        
        function res = eval(this, point1, point2)
           
            dist = sqrt((point1 - point2)' * (point1 - point2));
            res = this.Scale / dist;
            
        end
        
        function [res, num_evals] = eval_mat(this, data, rows, cols, filter)
           
            res = zeros(size(rows,2), size(cols,2));
                
            if (nargin == 4)
              
                if (numel(rows) < numel(cols)) 
                
                    for i = 1:numel(rows)

                        point_i = data(:,rows(i));
                        dists = sqrt(sum(bsxfun(@minus, data(:,cols), point_i).^2,1));
                        res(i,:) = 1.0 ./ dists;

                    end
                else
    
                    for i = 1:numel(cols)
                       
                        point_i = data(:,cols(i));
                        dists = sqrt(sum(bsxfun(@minus, data(:,rows), point_i).^2,1));
                        
                        res(:,i) = this.Scale ./ dists;
                        
                    end
                    
                end
                num_evals = numel(rows) * numel(cols);
               
            else
                % then we have a filter
                [row_inds,col_inds] = find(filter);

                for i = 1:numel(row_inds)

                    row_point = data(rows(row_inds(i)),:);
                    col_point = data(cols(col_inds(i)),:);

                    res(row_inds(i), col_inds(i)) = this.eval(row_point, col_point);

                end
                % this is the number of non zeros in the filter
                num_evals = numel(row_inds);
                
            end
        end
        
    end
    
    
    
end