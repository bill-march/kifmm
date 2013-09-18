classdef KIFMM

    
    properties
       
        Data
        
        Tree
        
        MaxTreeDepth
        
        MinNodeSize
        
        MaxLeafCount
        
        Kernel
        
        NumInterpolationPoints
        
        
    end
    
    
    methods
        
        % Main function for kernel-independent fmm
        function obj = kifmm(data, max_leaf_count, min_node_size, max_tree_depth, kernel, num_interpolation_points)

            obj.MaxLeafCount = max_leaf_count;
            obj.MinNodeSize = min_node_size;
            obj.MaxTreeDepth = max_tree_depth;
            
            obj.Tree = Octree1DClass(data, obj.MaxLeafCount, obj.MinNodeSize, obj.MaxTreeDepth);

            % make sure we get the sorted data
            obj.Data = obj.Tree.Data;
            
            obj.Kernel = kernel;
            
            obj.NumInterpolationPoints = num_interpolation_points;
            
            obj.PrecomputeRepresentation(obj);
            
        end

        function PrecomputeRepresentation(this)
 
            % First, go over all leaves
            
            num_leaves = size(this.Tree.LeafNodeList, 2);
            
            for i = 1:num_leaves
               
                leaf_ind = this.Tree.LeafNodeList(i);
                leaf_node = this.Tree.NodeList(leaf_ind);
                
                interpolation_points = this.Tree.SampleFarField(leaf_node, this.NumInterpolationPoints);
                num_interpolation_points = size(interpolation_points, 2);

                leaf_node.InterpolationPoints = interpolation_points;
                
                direct_matrix = zeros(this.NumInterpolationPoints, leaf_node.Count);
                direct_matrix_trans = zeros(this.NumInterpolationPoints, leaf_node.Count);

                for i = 1:num_interpolation_points

                   for j = leaf_node.Begin:leaf_node.End

                       direct_matrix(i,j) = kernel.eval(obj.Data(interpolation_points(i)), data(j));
                        % we're transposing this in place 
                        % note that for symmetric kernels, these matrices
                        % are each other's transpose
                       direct_matrix_trans(i,j) = kernel.eval(data(j), data(interpolation_points(i)));

                   end

               end

               [Acol, proj, skeleton] = InterpolativeDecomposition(direct_matrix);
               
               leaf_node.OutgoingSkeletonSize = size(skeleton, 2);
               leaf_node.OutgoingSkeleton = skeleton;
               leaf_node.ProjMatrix = proj;
               
               %also, compute the incoming representation
               [Brow_trans, eval_trans, iskel] = InterpolativeDecomposition(direct_matrix_trans);
               leaf_node.IncomingSkeletonSize = size(iskel, 2);
               leaf_node.IncomingSkeleton = iskel;
               leaf_node.EvalMatrix = eval_trans';
               
               
            end
            
            % now, go from finer to coarser levels
            % anything at the bottom level is always a leaf, so skip it
            for level = this.Tree.MaxDepth - 1:2
               
                num_nodes = size(this.Tree.NodesPerLevelList{level}, 2);
                
                for i = 1:num_nodes
                   
                    node_ind = this.Tree.NodesPerLevelList{level}(i);
                    node = this.Tree.NodeList(node_ind);
                    
                    % skip it if it is a leaf
                    if (~node.IsLeaf()) 
                        
                        % form the ii and oo operators
                        
                        left_ind = node.Index * 2;
                        right_ind = node.Index * 2 + 1;
                        
                        left_node = this.Tree.NodeList(left_ind);
                        right_node = this.Tree.NodeList(right_ind);
                        
                        % we can assume that they are both not empty, since
                        % leaves have already been taken care of 
                        if (left_node.is_empty()) 
                            
                            % my skeleton will just be the same as the 
                            % left_node?
                            % Problem: what if the interpolation points
                            % used to construct the child's representation
                            % are not valid for the parent?
                            % This isn't a problem, we will only be using
                            % the outgoing rep for things that are well
                            % separated from this node -- making it farther
                            % away won't make the ones we're using less
                            % valid
                            
                            node.CopySkeletons(right_node);
                            
                        elseif(right_node.is_empty())

                            node.CopySkeletons(left_node);
                            
                        else
                            
                            % need to merge them
                            
                            node.InterpolationPoints = this.Tree.SampleFarField(node, this.NumInterpolationPoints);
                            num_interpolation_points = size(node.InterpolationPoints, 2);
                            num_child_skeletons = left_node.OutgoingSkeletonSize + right_node.OutgoingSkeletonSize;

                            
                            A = zeros(num_interpolation_points, num_child_skeletons);
                            A_eval = zeros(num_interpolation_points, num_child_skeletons);
                            
                            for j = 1:num_interpolation_points
                               
                                j_point = this.Data(node.InterpolationPoints(j));
                                
                                for k = 1:left_node.OutgoingSkeletonSize
                                
                                    k_point = this.Data(left_node.OutgoingSkeleton(k));
                                    
                                    A(j,k) = this.Kernel.eval(j_point, k_point);
                                    A_eval(j,k) = this.Kernel.eval(k_point, j_point);
                                    
                                end

                                for k = 1:right_node.OutgoingSkeletonSize
                                
                                    k_point = this.Data(right_node.OutgoingSkeleton(k));
                                    
                                    A(j,k+left_node.OutgoingSkeletonSize) = this.Kernel.eval(j_point, k_point);
                                    A_eval(j,k+left_node.OutgoingSkeletonSize) = this.Kernel.eval(k_point, j_point);
                                    
                                end

                            end
                            
                            % now, we've formed the matrix A, so decompose
                            % it
                            [this_Acol, Z, skel] = InterpolativeDecomposition(A);
                            node.ProjMatrix = Z;
                            node.OutgoingSkeleton = zeros(1, sizeof(skel,2));
                            
                            [this_Brow, Zeval, iskel] = InterpolativeDecomposition(A_eval);
                            node.EvalMatrix = Zeval';
                            node.IncomingSkeleton = zeros(1, sizeof(iskel,2));
                            
                            
                            for j = 1:size(skel,2)
                               
                                this_ind = skel(j);
                                if (this_ind <= left_node.OutgoingSkeletonSize)
                                    node.OutgoingSkeleton(j) = left_node.OutgoingSkeleton(this_ind);
                                else
                                    node.OutgoingSkeleton(j) = right_node.OutgoingSkeleton(this_ind);                                    
                                end
                                
                            end
                            
                            for j = 1:size(iskel,2)
                               
                                this_ind = iskel(j);
                                if (this_ind <= left_node.IncomingSkeletonSize)
                                    node.IncomingSkeleton(j) = left_node.IncomingSkeleton(this_ind);
                                else
                                    node.IncomingSkeleton(j) = right_node.IncomingSkeleton(this_ind);                                    
                                end
                                
                            end
                            
                        end   
                        
                    end
                    
                end
                
            end
            
            % now, form the oi operators
            
            num_nodes = size(this.Tree.NodeList, 2);
            
            for i = 1:num_nodes
               
                this_node = this.Tree.NodeList(i);
                
                num_interacting = size(this_node.InteractionList, 2);
                
                this_node.OIMatrices = cell([1, num_interacting]);
                
                for j = 1:num_interacting
                   
                    int_ind = this_node.InteractionList(j);
                    int_node = this.Tree.NodeList(int_ind);
                    
                    % now, compute the oi operator between this_node and
                    % int_node
                    
                    this_oi = zeros(this_node.IncomingSkeletonSize, int_node.OutgoingSkeletonSize);
                    
                    for k = 1:this_node.IncomingSkeletonSize
                       
                        k_point = this.Data(this_node.IncomingSkeleton(k));
                        
                        for l = 1:int_node.OutgoingSkeletonSize
                        
                            l_point = this.Data(int_node.OutgoingSkeleton(l));
                            
                            this_oi(k,l) = this.Kernel.eval(k_point, l_point);
                            
                        end
                        
                    end
                    
                    this_node.OIMatrices{j} = this_oi;
                    
                end
                
            end
            
               
       
        end

    end

end




