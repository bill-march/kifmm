classdef Tree

    properties
        NumLevels % number of levels in the tree
        Nodes % master list of nodes, 0 is the root
        NodesPerLevel % the list of nodes in each level, these need to be references
        
    end
    
    methods
       
        function obj=Tree(Data, MinBoxSize, NumLevels)
            RootBox = ComputeBounds(Data);
            root = TreeNode(bounds, 1, Data.ncols(), NULL);
            nodes.push_back(root);
            nodes_per_level[1].push_back(root);
            
        end
        
        function [new_nodes] = SplitNode(this, node_id)
           
            if nodes[node_id].Count > MinNodeSize
                
                new_nodes = PartitionMatrix(Data, nodes[node_id].begin, nodes[node_id].count);
                
            end % do we have enough nodes to fool with
            
        end
            
    end

end

function [new_nodes] = PartitionMatrix(data, begin, count)

    PivotVal = Bounds.midpoint;
    
    for d=1:dimension
    
        for i=begin:begin+count

            point = data[i];
            

        end

    end


end


        




classdef TreeNode
   
    properties
       Bounds % the bounding box 
       Parent % id of the parent node
       Children % list of ids of child nodes
       Begin % index of the first point
       Count % number of points
       OutgoingRep
       IncomingRep
       
    end
    
    methods
       
        function obj=TreeNode(bounds, begin, count, parent)
            Bounds = bounds;
            Begin = begin;
            Count = count;
            Parent = parent;
        end
        
        
        
    end
    
    
end

classdef BoundingBox
   
    properties
        MinValues
        MaxValues
    end
    
    
    
end
