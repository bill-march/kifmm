classdef Octree1DClass < handle
    
    properties
        
        % children of node i are at 2i, 2i+1
        NodeList % contains all the nodes
        NodesPerLevelList % contains all nodes at each level
        
        % the total depth of the tree
        Depth
        
        % the inputs
        Data
        
        % the list of all leaf nodes
        LeafNodeList
        
        % the size of the box bounding the root node
        NodeSize
        
        % paramters for tree building
        MaxLeafSize
        MinLeafWidth
        MaxDepth
        
    end
    
    methods
        
        function obj = Octree1DClass(Data, MaxLeafSize, MinLeafWidth, MaxDepth)
            
            obj.Depth = 1;
            
            num_points = size(Data, 2);
            
            obj.MaxLeafSize = MaxLeafSize;
            obj.MinLeafWidth = MinLeafWidth;
            obj.MaxDepth = MaxDepth;
            
            % Preallocate the list of nodes 
            %obj.NodeList(2^MaxDepth) = OctreeNode;
            obj.NodeList = OctreeNode.empty(2^MaxDepth,0);
            obj.LeafNodeList = OctreeNode.empty(ceil(num_points / MaxLeafSize), 0);
            
            obj.Data = sort(Data);
            
            root_node = OctreeNode(1, num_points, 1, 1, obj.Data(1), obj.Data(end));
            
            obj.NodeSize = obj.Data(end) - obj.Data(1);
            
            SplitNode(obj, root_node);
            
            FillInLists(obj);
            
            
        end
        
        
        function SplitNode(this, node)
            
            this.NodeList(node.Index) = node;
            this.NodesPerLevelList{node.Depth}(end+1) = node;
            
            % fill me in
            node_width = node.MaxVal - node.MinVal;

            % can we split the node at all?
            if ((node.Count > this.MaxLeafSize) && (node_width > this.MinLeafWidth) && (node.Depth < this.MaxDepth)) 
                
                split_val = (node.MinVal + node.MaxVal) / 2;

                % should be the index of the last point in the left child
                % TODO: implement this non-stupidly
                begin = node.Begin;
                endval = node.End;
                split_index = find((this.Data(:, node.Begin:node.End) < split_val), 1, 'last');
                % if the split index is empty, then all the points go on
                % the right side
                if (isempty(split_index)) 
                    split_index = 0;
                end
                
                left_index = 2 * node.Index;
                right_index = 2 * node.Index + 1;
                
                left_endpoint = node.MinVal + 2^(-node.Depth) * this.NodeSize;
                
                left_count = split_index;
                right_begin = node.Begin + left_count;
                right_count = node.Count - split_index;
                
                
                left_child = OctreeNode(node.Begin, left_count, node.Depth+1, left_index, node.MinVal, left_endpoint);
                right_child = OctreeNode(right_begin, right_count, node.Depth+1, right_index, left_endpoint, node.MaxVal);
                
                this.NodeList(left_index) = left_child;
                this.NodeList(right_index) = right_child;
                
                % only try to split them if they have any points
                if (left_count > 0)
                    node.HasLeft = true;
                    SplitNode(this, left_child);
                end
                
                if (right_count > 0) 
                    node.HasRight = true;
                    SplitNode(this, right_child);
                end
                
            else 
                % this node is a leaf
                this.LeafNodeList(end+1) = node;
                this.Depth = max([this.Depth node.Depth]);
                
            end
            
        end
        
        
        function FillInLists(this, node)
        
            if (node.IsLeaf())
               
                
                
            end
            
            if (node.HasLeft) 
                left_child = this.NodeList(2 * node.Index);
                this.FillInLists(left_child);
            end
            
            if (node.HasRight)
                right_child = this.NodeList(2 * node.Index + 1);
                this.FillInLists(right_child);
            end  
            
        end
        
        
        function Print(this, node)
           
            node.Print();
            
            if (node.HasLeft) 
                left_child = this.NodeList(2 * node.Index);
                this.Print(left_child);
            end
            
            if (node.HasRight)
                right_child = this.NodeList(2 * node.Index + 1);
                this.Print(right_child);
            end
            
        end
        
    end
    
end

