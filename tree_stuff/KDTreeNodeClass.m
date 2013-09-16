classdef KDTreeNodeClass
    
    properties
    
        Bounds % the bounding box
        Children % indices of the children
        Begin % the data index where this node begins
        Count % the number of data points in this node
        Last %
        
    end
    
    methods
        
        function obj = KDTreeNodeClass(begin, count, bounds)
        
            obj.Bounds = bounds;
            obj.Begin = begin;
            obj.Count = count;
            obj.Last = begin + count;
            
        end
        
    end
    
end

