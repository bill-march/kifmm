classdef KIFMM < handle

    
    properties
       
        Data
        NumPoints
        
        Charges
        Potentials
        
        Tree % the tree class
        
        MaxTreeDepth
        
        MinNodeSize
        
        MaxLeafCount
        
        Kernel
        
        
    end
    
    
    methods
        
        % Main function for kernel-independent fmm
        function obj = KIFMM(data, max_leaf_count, min_node_size, max_tree_depth, kernel)

            if (nargin > 0)
                obj.MaxLeafCount = max_leaf_count;
                obj.MinNodeSize = min_node_size;
                obj.MaxTreeDepth = max_tree_depth;

                obj.Tree = Octree1DClass(data, max_leaf_count, min_node_size, max_tree_depth);

                % make sure we get the sorted data
                obj.Data = obj.Tree.Data;
                obj.NumPoints = size(obj.Data, 2);

                obj.Kernel = kernel;

                obj.PrecomputeRepresentation();
            end
        end

        function potentials = ComputePotentials(this, charges)
            
            % First, permute the charges in the same way we permuted the
            % points in the tree
            this.Charges = zeros(this.NumPoints,1);
            this.Potentials = zeros(this.NumPoints,1);
            
            for i = 1:this.NumPoints
                this.Charges(i) = charges(this.Tree.NewFromOld(i));
            end
            
            % for all leaves, compute the outgoing rep
            num_leaves = size(this.Tree.LeafNodeList, 2);
            
            for i = 1:num_leaves
               
                leaf_ind = this.Tree.LeafNodeList(i);
                leaf_node = this.Tree.NodeList(leaf_ind);
                
                leaf_node.PsiVector = leaf_node.ProjMatrix * this.Charges(leaf_node.Begin:leaf_node.End);
                
            end
            
            % Upward pass
            
            for level = this.Tree.MaxDepth - 1:-1:3
               
                num_nodes = size(this.Tree.NodesPerLevelList{level}, 2);
                
                for i = 1:num_nodes
                   
                    node_ind = this.Tree.NodesPerLevelList{level}(i);
                    node = this.Tree.NodeList(node_ind);
                    
                    if (~node.IsLeaf()) 
                        
                        left_ind = node.Index * 2;
                        right_ind = node.Index * 2 + 1;
                        
                        left_node = this.Tree.NodeList(left_ind);
                        right_node = this.Tree.NodeList(right_ind);
                        
                        node.PsiVector = zeros(node.OutgoingSkeletonSize,1);
                        
                        if (~left_node.is_empty()) 
                            %node.PsiVector = node.PsiVector + node.ProjMatrix(:,1:node.NumOutLeft) * this.Data(left_node.OutgoingSkeleton) * this.Charges(left_node.OutgoingSkeleton);
                            node.PsiVector = node.PsiVector + node.ProjMatrix(:,1:node.NumOutLeft) * left_node.PsiVector;
                        end
                        
                        if (~right_node.is_empty())
                            %node.PsiVector = node.PsiVector + node.ProjMatrix(:, node.NumOutLeft+1:end) * this.Data(right_node.OutgoingSkeleton) * this.Charges(right_node.OutgoingSkeleton);
                            node.PsiVector = node.PsiVector + node.ProjMatrix(:,node.NumOutLeft+1:end) * right_node.PsiVector;
                        end
                        
                    end
                    
                end
                
            end % end of upward pass
            
            % Downward pass
            
            this.Tree.NodeList(1).PhiVector = zeros(this.NumPoints, 1);
            this.Tree.NodeList(2).PhiVector = zeros(this.Tree.NodeList(2).Count, 1);
            this.Tree.NodeList(3).PhiVector = zeros(this.Tree.NodeList(3).Count, 1);
            
            this.Tree.NodeList(2).EvalMatrix = eye(this.Tree.NodeList(2).Count);
            this.Tree.NodeList(3).EvalMatrix = eye(this.Tree.NodeList(3).Count);
            
            for level = 3:this.Tree.MaxDepth
               
                num_nodes = size(this.Tree.NodesPerLevelList{level},2);
                
                for i = 1:num_nodes
                   
                    node_ind = this.Tree.NodesPerLevelList{level}(i);
                    node = this.Tree.NodeList(node_ind);
                    
                    parent_ind = floor(node.Index / 2);
                    parent_node = this.Tree.NodeList(parent_ind);

                    % am I the left or right child of my parent?
                    is_left = (mod(node.Index, 2) == 0);
                    
                    if (is_left) 
                        node.PhiVector = parent_node.EvalMatrix(1:node.Count,:) * parent_node.PhiVector;
                    else
                        node.PhiVector = parent_node.EvalMatrix(node.Begin - parent_node.Begin + 1:end, :) * parent_node.PhiVector;
                    end
                    
                    num_interacting = size(node.InteractionList, 2);
                    
                    for j = 1:num_interacting
                       
                        int_ind = node.InteractionList(j);
                        int_node = this.Tree.NodeList(int_ind);
                        
                        oi_mat = node.OIMatrices{j};
                        
                        node.PhiVector = node.PhiVector + oi_mat * int_node.PsiVector;
                        
                    end
                    
                end
                
            end % end of downward pass
            
            % now, finish up at the leaves
            
            num_leaves = size(this.Tree.LeafNodeList, 2);
            
            for i = 1:num_leaves
            
                leaf_ind = this.Tree.LeafNodeList(i);
                leaf_node = this.Tree.NodeList(leaf_ind);
                
                this.Potentials(leaf_node.Begin:leaf_node.End) = leaf_node.EvalMatrix * leaf_node.PhiVector;
                
                for point_ind = leaf_node.Begin:leaf_node.End
                   
                    point = this.Data(point_ind);
                    
                    num_neighbors = size(leaf_node.NearFieldList, 2);
                    for n_node_ind = 1:num_neighbors
                       
                        n_node = this.Tree.NodeList(leaf_node.NearFieldList(n_node_ind));
                        
                        for n_point_ind = n_node.Begin:n_node.End
                           
                            n_point = this.Data(n_point_ind);
                            n_charge = this.Charges(n_point_ind);
                            
                            this.Potentials(point_ind) = this.Potentials(point_ind) + this.Kernel.eval(point, n_point) * n_charge;
                            
                        end
                        
                    end
                    
                    % now, do the self-interactions
                    for op_ind = leaf_node.Begin:leaf_node.End
                       
                        op = this.Data(op_ind);
                        op_charge = this.Charges(op_ind);
                        this.Potentials(point_ind) = this.Potentials(point_ind) + this.Kernel.eval(point, op) * op_charge;
                        
                    end
                    
                end
                
            end
            
            % now, unpermute the potentials for output
            potentials = zeros(this.NumPoints, 1);
            
            this.Tree.NewFromOld
            
            for i = 1:this.NumPoints
               
                potentials(this.Tree.NewFromOld(i)) = this.Potentials(i);
                
            end
            
        end
        
        
        function PrecomputeRepresentation(this)
 
            % First, go over all leaves
            
            num_leaves = size(this.Tree.LeafNodeList, 2);
            
            for i = 1:num_leaves
               
                leaf_ind = this.Tree.LeafNodeList(i);
                leaf_node = this.Tree.NodeList(leaf_ind);
                
                interpolation_points = this.Tree.SampleFarField(leaf_node);

               [direct_matrix, proj, skeleton] = InterpolativeDecomposition(this.Data, interpolation_points, leaf_node.Begin:leaf_node.End, this.Kernel);
               
               leaf_node.OutgoingSkeletonSize = size(skeleton, 2);
               leaf_node.OutgoingSkeleton = skeleton;
               leaf_node.ProjMatrix = proj;
               
               %also, compute the incoming representation
               % IMPORTANT: this doesn't allow for non-symmetric kernels 
               [direct_matrix_trans, eval_trans, iskel] = InterpolativeDecomposition(this.Data, interpolation_points, leaf_node.Begin:leaf_node.End, this.Kernel);
               leaf_node.IncomingSkeletonSize = size(iskel, 2);
               leaf_node.IncomingSkeleton = iskel;
               leaf_node.EvalMatrix = eval_trans';
               
               
            end
            
            % now, go from finer to coarser levels
            % anything at the bottom level is always a leaf, so skip it
            for level = this.Tree.MaxDepth - 1:-1:3
               
                num_nodes = size(this.Tree.NodesPerLevelList{level}, 2);
                
                for i = 1:num_nodes
                   
                    node_ind = this.Tree.NodesPerLevelList{level}(i);
                    node = this.Tree.NodeList(node_ind);
                    
                    % skip it if it is a leaf, we already formed the
                    % matrices we need
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
                            %
                            % nothing owned by the parent can be well
                            % separated from the child
                            
                            node.CopySkeletons(right_node);
                            node.NumOutRight = right_node.OutgoingSkeletonSize;
                            node.NumInRight = right_node.IncomingSkeletonSize;
                            
                        elseif(right_node.is_empty())

                            node.CopySkeletons(left_node);
                            node.NumOutLeft = left_node.OutgoingSkeletonSize;
                            node.NumInLeft = left_node.IncomingSkeletonSize;
                            
                        else
                            
                            % need to merge them
                            
                            interpolation_points = this.Tree.SampleFarField(node);
                            
                            % we're assuming that the skeletons are
                            % disjoint because the nodes are disjoint
                            %merged_skeletons = union(left_node.OutgoingSkeleton, right_node.OutgoingSkeleton);
                            merged_skeletons = [left_node.OutgoingSkeleton, right_node.OutgoingSkeleton];
                            
                            merged_out_skel = [left_node.IncomingSkeleton, right_node.IncomingSkeleton];
                            
                            % now, we've formed the matrix A, so decompose
                            % it
                            [this_Acol, Z, skel] = InterpolativeDecomposition(this.Data, interpolation_points, merged_skeletons, this.Kernel);
                            node.ProjMatrix = Z;
                            node.OutgoingSkeleton = skel;
                            node.OutgoingSkeletonSize = size(skel,2);
                            
                            [this_Brow, Zeval, iskel] = InterpolativeDecomposition(this.Data, interpolation_points, merged_out_skel, this.Kernel);
                            node.EvalMatrix = Zeval';
                            node.IncomingSkeleton = iskel;
                            node.IncomingSkeletonSize = size(iskel,2);
                            
                            % need to set the indices in the node
                            [row, first_right_ind_out] = find(node.OutgoingSkeleton > left_node.End, 1, 'first');
                            [row, first_right_ind_in] = find(node.IncomingSkeleton > left_node.End, 1, 'first');
                            
                            node.NumOutLeft = first_right_ind_out - 1;
                            node.NumOutRight = node.OutgoingSkeletonSize - node.NumOutLeft;
                            
                            node.NumInLeft = first_right_ind_in - 1;
                            node.NumInRight = node.IncomingSkeletonSize - node.NumInLeft;
                            
                        end   
                        
                    end
                    
                end
                
            end
            
            % now, form the oi operators
            
            num_nodes = size(this.Tree.NodeList, 2);
            
            for i = 1:num_nodes
               
                this_node = this.Tree.NodeList(i);
                
                if (~this_node.is_empty()) 

                    num_interacting = size(this_node.InteractionList, 2);
                
                    %this_node.OIMatrices = cell([1, num_interacting]);

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
        
        function potentials = ComputePotentialsNaive(this, data, charges)
           
            num_points = size(data, 2);
            potentials = zeros(num_points, 1);
            
            for i = 1:num_points
                
                this_potential = 0;
                
                for j = 1:num_points
                   
                    this_ker = this.Kernel.eval(data(i), data(j));
                    this_potential = this_potential + this_ker * charges(j);
                    
                end
                
                potentials(i) = this_potential;
                
            end
            
        end

    end

end




