classdef UniformSampler
    % Samples from the far-field uniformly at random (iid)
    
    properties
        
        Probability
        
    end
    
    methods

        function obj = UniformSampler(prob)
    
            obj.Probability = prob;
            
        end
        
        function filter = GenFilter(this, num_rows, num_cols)
           
            filter = zeros(num_rows, num_cols);
            
            for i = 1:num_rows
               
                for j = 1:num_cols
                   
                    p = rand(1);
                    
                    if (p < this.Probability) 
                        filter(i,j) = 1;
                    end
                    
                end
                
            end
            
        end
                
        function [filter] = SampleFarField(this, row_inds, col_inds)
            
            num_rows = size(row_inds,2);
            num_cols = size(col_inds,2);

            % this doesn't generate a dense matrix when the probability is
            % 1
            %filter = sprand(num_rows, num_cols, this.Probability);
            
            %num_rows
            %num_cols
            %this.Probability
            filter = this.GenFilter(num_rows, num_cols);
            %nnz(filter)
            
        end
        
    end
    
end

