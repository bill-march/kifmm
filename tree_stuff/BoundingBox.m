classdef BoundingBox
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        LowerBounds
        UpperBounds
        
    end
    
    methods
        
        function obj = BoundingBox(data)
        
            % should give me the min in each row over all columns, which is
            % what I want
            obj.LowerBounds = min(data,[],2)
            obj.UpperBounds = max(data,[],2)
            
        end

        function largest = SplitDim(this)
           
            [max_val, argmax] = max(this.UpperBounds - this.LowerBounds);
            largest = argmax;
            
        end
        
        function dist = MinDistanceSq(this, other)
        
            lower = this.LowerBounds - other.UpperBounds;
            upper = other.LowerBounds - this.UpperBounds;
            
            % using the trick in mlpack
            vec = lower + abs(lower) + upper + abs(upper);
            
            % now, take the dot product
            dist = vec' * vec;
            
        end
        
        function [width, dim] = MaxWidth(this)
           
            [width,dim] = max(this.UpperBounds - this.LowerBounds);
            
        end
        
        
        
    end
    
end

