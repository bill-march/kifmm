classdef OctreeNode < handle
    
    properties
        IncomingRepresentation
        OutgoingRepresentation
    
        % the index in the global array of the first point owned by the
        % node
        Begin

        % the number of points owned by the node
        Count
        
        % the index in the global data matrix of the last point owned by 
        % the node
        End
        
        HasLeft
        HasRight
        
        MinVal
        MaxVal
        
        Depth
        
        % this node's index in the global list of nodes
        Index
        
        % this nodes lists 
        InteractionList
        NearFieldList
    
    end
    
    methods
        
        function obj = OctreeNode(begin, count, depth, index, min_val, max_val)
            
            if (nargin > 0)
                if (count > 0) 
                    obj.Begin = begin;
                    obj.Count = count;
                    obj.Depth = depth;
                    obj.Index = index;
                    obj.MinVal = min_val;
                    obj.MaxVal = max_val;
                    obj.End = begin + count - 1;
                    obj.HasLeft = false;
                    obj.HasRight = false;
                end
            end            
        end
        
        function Print(this)
           
            'Node:'
            begin = this.Begin
            count = this.Count
            
        end
                
    end
    
end

