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
        
        SkeletonSize
        
        
    end
    
    
    methods
        
        % Main function for kernel-independent fmm
        function obj = KIFMM(data, max_leaf_count, min_node_size, max_tree_depth, skeleton_size, kernel)

            if (nargin > 0)
                obj.MaxLeafCount = max_leaf_count;
                obj.MinNodeSize = min_node_size;
                obj.MaxTreeDepth = max_tree_depth;
                
                obj.SkeletonSize = skeleton_size;

                obj.Tree = Octree1DClass(data, max_leaf_count, min_node_size, max_tree_depth);

                % make sure we get the sorted data
                obj.Data = obj.Tree.Data;
                obj.NumPoints = size(obj.Data, 2);

                obj.Kernel = kernel;

                obj.PrecomputeRepresentation();
            end
        end
        
         
        function [proj, skeleton] = DecomposeKernel(this, row_inds, col_inds)

            Acol = this.Kernel.eval_mat(this.Data, row_inds, col_inds);
            
            [proj, skeleton] = InterpolativeDecomposition(Acol, this.SkeletonSize);
            
        end
        
        function [eval, skeleton] = DecomposeKernelTrans(this, row_inds, col_inds)
           
            Btrans = this.Kernel.eval_mat(this.Data, col_inds, row_inds);
            [proj,skeleton] = InterpolativeDecomposition(Btrans, this.SkeletonSize);
            eval = proj';
            
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
                            node.PsiVector = node.PsiVector + node.ProjMatrix(:,1:left_node.OutgoingSkeletonSize) * left_node.PsiVector;
                        end
                        
                        if (~right_node.is_empty())
                            node.PsiVector = node.PsiVector + node.ProjMatrix(:,end-right_node.OutgoingSkeletonSize+1:end) * right_node.PsiVector;
                        end
                        
                    end
                    
                end
                
            end % end of upward pass
            
            % Downward pass
            
            % initialize the top of the tree to be empty
            this.Tree.NodeList(1).PhiVector = zeros(this.SkeletonSize,1);
            
            this.Tree.NodeList(2).PhiVector = zeros(this.SkeletonSize, 1);
            this.Tree.NodeList(3).PhiVector = zeros(this.SkeletonSize, 1);
            
            this.Tree.NodeList(2).EvalMatrix = eye(2 * this.SkeletonSize, this.SkeletonSize);
            this.Tree.NodeList(3).EvalMatrix = eye(2 * this.SkeletonSize, this.SkeletonSize);
            
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
                        node.PhiVector = parent_node.EvalMatrix(1:node.IncomingSkeletonSize,:) * parent_node.PhiVector;
                    else
                        node.PhiVector = parent_node.EvalMatrix(end - node.IncomingSkeletonSize + 1:end,:) * parent_node.PhiVector;
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
            
            potentials(this.Tree.NewFromOld) = this.Potentials;
            
        end
        
        
        function PrecomputeRepresentation(this)
 
            % First, go over all leaves
            
            num_leaves = size(this.Tree.LeafNodeList, 2);
            
            for i = 1:num_leaves
               
                leaf_ind = this.Tree.LeafNodeList(i);
                leaf_node = this.Tree.NodeList(leaf_ind);
                
                interpolation_points = this.Tree.SampleFarField(leaf_node);

               [proj, skeleton] = this.DecomposeKernel(interpolation_points, leaf_node.Begin:leaf_node.End);
               
               leaf_node.OutgoingSkeletonSize = size(skeleton, 2);
               [leaf_node.OutgoingSkeleton, sort_id] = sort(skeleton + leaf_node.Begin - 1);
               leaf_node.ProjMatrix = proj(sort_id,:);
               
               %also, compute the incoming representation
               [eval, iskel] = this.DecomposeKernelTrans(leaf_node.Begin:leaf_node.End, interpolation_points);
               leaf_node.IncomingSkeletonSize = size(iskel, 2);
               [leaf_node.IncomingSkeleton, sort_id] = sort(iskel + leaf_node.Begin - 1);
               leaf_node.EvalMatrix = eval(:,sort_id);
               
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
                    % this also checks that it's not empty
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
                            
                            merged_in_skel = [left_node.IncomingSkeleton, right_node.IncomingSkeleton];
                            
                            % now, we've formed the matrix A, so decompose
                            % it
                            [Z, skel] = this.DecomposeKernel(interpolation_points, merged_skeletons);
                            [node.OutgoingSkeleton, sort_id] = sort(merged_skeletons(skel));
                            node.ProjMatrix = Z(sort_id,:);
                            node.OutgoingSkeletonSize = size(skel,2);
                            
                            [Zeval, iskel] = this.DecomposeKernelTrans(merged_in_skel, interpolation_points);
                            [node.IncomingSkeleton, sort_id] = sort(merged_in_skel(iskel));
                            node.EvalMatrix = Zeval(:,sort_id);
                            node.IncomingSkeletonSize = size(iskel,2);
                            
                            % need to set the indices in the node
                            [~, first_right_ind_out] = find(node.OutgoingSkeleton > left_node.End, 1, 'first');
                            [~, first_right_ind_in] = find(node.IncomingSkeleton > left_node.End, 1, 'first');
                            
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
                
                    for j = 1:num_interacting

                        int_ind = this_node.InteractionList(j);
                        int_node = this.Tree.NodeList(int_ind);

                        % now, compute the oi operator between this_node and
                        % int_node
                        this_oi = this.Kernel.eval_mat(this.Data, this_node.IncomingSkeleton, int_node.OutgoingSkeleton);
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




