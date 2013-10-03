classdef GaussianKernel < handle
   
    properties
        
        Bandwidth
        Scale 
        
    end
    
    
    methods
       
        function obj = GaussianKernel(band, scale)
            
            obj.Bandwidth = band;
            
            if (nargin > 1)
                obj.Scale = scale;
            else 
                obj.Scale = 1;
            end
            
        end
        
        function res = eval(this, point1, point2)
           
            dist_sqr = (point1 - point2)^2;
            res = exp(-1 * dist_sqr * this.Bandwidth);
            
        end
        
        function [res, num_evals] = eval_mat(this, data, rows, cols, filter)
           
            res = zeros(size(rows,2), size(cols,2));
                
            if (nargin == 4)
              
                for i = 1:numel(rows)
                   
                    point_i = data(rows(i));
                    
                    dists = (data(cols) - point_i).^2;
                    res(i,:) = exp(-dists * this.Bandwidth);
                    
                end
                
%                 for i = 1:numel(rows)
%                
%                     point_i = data(rows(i));
% 
%                     for j = 1:numel(cols)
% 
%                         point_j = data(cols(j));
% 
%                         res(i,j) = this.eval(point_i,point_j);
%                         
%                     end
% 
%                 end
%                 
%                 num_evals = numel(rows) * numel(cols);
                
                
                
            
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