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
        
        % the permutation used in sorting
        NewFromOld
        
    end
    
    methods
        
        function obj = Octree1DClass(Data, MaxLeafSize, MinLeafWidth, MaxDepth)
            
            if (nargin > 0)
            
                obj.Depth = 1;

                num_points = size(Data, 2);

                obj.MaxLeafSize = MaxLeafSize;
                obj.MinLeafWidth = MinLeafWidth;
                obj.MaxDepth = MaxDepth;

                % Preallocate the list of nodes 
                %obj.NodeList(2^MaxDepth) = OctreeNode;
                obj.NodeList = OctreeNode.empty(2^MaxDepth,0);
                obj.LeafNodeList = int32.empty(ceil(num_points / MaxLeafSize), 0);
                obj.NodesPerLevelList = cell([1, MaxDepth]);

                [obj.Data, obj.NewFromOld] = sort(Data);

                root_node = OctreeNode(1, num_points, 1, obj.Data(1), obj.Data(end), 1);

                obj.NodeSize = obj.Data(end) - obj.Data(1);

                SplitNode(obj, root_node);

                FillInLists(obj, root_node, root_node);

            end
            
        end
        
        
        function SplitNode(this, node)
            
            this.NodeList(node.Index) = node;
            this.NodesPerLevelList{node.Depth}(end+1) = node.Index;
            
            % fill me in
            node_width = node.MaxVal - node.MinVal;

            % can we split the node at all?
            if ((node.Count > this.MaxLeafSize) && (node_width > this.MinLeafWidth) && (node.Depth < this.MaxDepth)) 
                
                split_val = (node.MinVal + node.MaxVal) / 2;

                split_index = find((this.Data(:, node.Begin:node.End) < split_val), 1, 'last');
                % if the split index is empty, then all the points go on
                % the right side
                if (isempty(split_index)) 
                    split_index = 0;
                end
                
                left_endpoint = node.MinVal + 2^(-node.Depth) * this.NodeSize;
                
                left_count = split_index;
                right_begin = node.Begin + left_count;
                right_count = node.Count - split_index;
                
                left_index = 2 * node.Index;
                right_index = left_index + 1;
                
                left_child = OctreeNode(node.Begin, left_count, node.Depth+1, node.MinVal, left_endpoint, left_index);
                right_child = OctreeNode(right_begin, right_count, node.Depth+1, left_endpoint, node.MaxVal, right_index);
                
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
                this.LeafNodeList(end+1) = node.Index;
                this.Depth = max([this.Depth node.Depth]);
                
            end
            
        end
        
        function FillInLists(this, node1, node2)
           
            if ((node1.Depth == node2.Depth) && node1.WellSeparated(node2))
                
                node1.InteractionList(end+1) = node2.Index;
            
                % don't add a self interaction here
            elseif (node1.IsLeaf() && node2.IsLeaf() && node1.Index ~= node2.Index) 
            
                node1.NearFieldList(end+1) = node2.Index;
                
            elseif(node1.IsLeaf()) 
                
                if (node2.HasLeft) 
                    left_child2 = node2.Index * 2;
                    this.FillInLists(node1, this.NodeList(left_child2));
                end

                if (node2.HasRight)
                    right_child2 = node2.Index * 2 + 1;
                    this.FillInLists(node1, this.NodeList(right_child2));
                end

            elseif(node2.IsLeaf())
                
                if (node1.HasLeft) 
                    left_child1 = node1.Index * 2;
                    this.FillInLists(this.NodeList(left_child1), node2);
                end

                if (node1.HasRight)
                    right_child1 = node1.Index * 2 + 1;
                    this.FillInLists(this.NodeList(right_child1), node2);
                end
                
            else % both are internal nodes
                
                if (node1.HasLeft) 
                    left_child1 = node1.Index * 2;
                    if (node2.HasLeft) 
                        left_child2 = node2.Index * 2;
                        this.FillInLists(this.NodeList(left_child1), this.NodeList(left_child2));
                    end

                    if (node2.HasRight)
                        right_child2 = node2.Index * 2 + 1;
                        this.FillInLists(this.NodeList(left_child1), this.NodeList(right_child2));
                    end
                end

                if (node1.HasRight) 
                    right_child1 = node1.Index * 2 + 1;
                    if (node2.HasLeft) 
                        left_child2 = node2.Index * 2;
                        this.FillInLists(this.NodeList(right_child1), this.NodeList(left_child2));
                    end

                    if (node2.HasRight)
                        right_child2 = node2.Index * 2 + 1;
                        this.FillInLists(this.NodeList(right_child1), this.NodeList(right_child2));
                    end
                end
                
            end
            
        end
                

        
        function Print(this)
           
            for depth = 1:this.MaxDepth
               
                num_nodes = size(this.NodesPerLevelList{depth}, 2);
                
                for node_ind = 1: num_nodes
                   
                    cur_node_ind = this.NodesPerLevelList{depth}(node_ind);
                    cur_node = this.NodeList(cur_node_ind)
                    
                end
                
            end

        end
        
        % this should mostly work
        function [InterpolationInds] = SampleFarField(this, node, num_points)
            
            num_sampled = 0;
            
            width = node.MaxVal - node.MinVal;
            
            total_points = size(this.Data, 2);
            InterpolationInds = zeros(1,num_points);
            possible_inds = 1:total_points;
            num_possible_points = total_points;
            
            while (num_sampled < num_points)
            
                sample_ind = randi(num_possible_points);
                sample_point = this.Data(possible_inds(sample_ind));
                dist = node.MinDistance(sample_point);
                
                if (dist >= width) 
                    InterpolationInds(num_sampled+1) = possible_inds(sample_ind);
                    num_sampled = num_sampled + 1;
                end
                
                temp = possible_inds(total_points - num_sampled);
                possible_inds(total_points - num_sampled) = possible_inds(sample_ind);
                possible_inds(sample_ind) = temp;
                
                num_possible_points = num_possible_points - 1;
                
            end
            
        end
        
    end
    
end

