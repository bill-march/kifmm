classdef UniformSampler
    % Samples from the far-field uniformly at random (iid)
    % this is based on the method of Achlioptas et al.
    
    properties
        
        Probability
        
    end
    
    methods

        function obj = UniformSampler(prob)
    
            obj.Probability = prob;
            
        end
        
        function filter = GenFilter(this, num_rows, num_cols)
           
            p = rand(num_rows, num_cols);
            filter = ~ceil(p - this.Probability);
            
        end
                
        function [filter] = SampleFarField(this, row_inds, col_inds)
            
            num_rows = numel(row_inds);
            num_cols = numel(col_inds);

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

