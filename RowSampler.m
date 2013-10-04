classdef RowSampler
    % Samples rows independently and uniformly at random
    
    properties
        
        Probability
        
    end
    
    methods

        function obj = RowSampler(prob)
    
            obj.Probability = prob;
            
        end
   
        function rows_out = GenFilter(this, rows_in)
           
            p = rand(1,numel(rows_in));
            filter = ~ceil(p - this.Probability);
            
            rows_out = rows_in(filter > 0);
            
        end
        
            
    end
    
end

