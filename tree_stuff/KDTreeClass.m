classdef KDTreeClass
   

    properties
        
        % list of nodes, root at 0
        Nodes
        % the data points (permuted)
        Data
        % The permutation
        OldFromNew
        
    end
    
    
    methods
        
        % constructor
        function obj = KDTreeClass(Data, MaxLeafSize, MinNodeExtent) 
        
            Nodes(1) = KDTreeClass;
            num_data = size(Data, 2);
            obj.Data = Data;
            
            box = BoundingBox(Data);
            
            root_node = KDTreeNode(1,num_data,box);
            Nodes(1) = root_node;
            
            SplitNode(root_node);
            
        end
        
        
        function SplitNode(this, node, MaxLeafSize, MinNodeExtent)
        
            % check that we can split at all
            if (node.Count > MaxLeafSize && node.Bounds.MaxWidth > MinNodeExtent) 
                
                splitdim = node.Bounds.MaxWidth();
                
                % IMPORTANT: this might be copying, which I don't want
                node_data = this.Data[~,node.Begin:node.Last];
                
                
                
                
                
            end
            
        end
        
            
            
    end
    
    
    
end

