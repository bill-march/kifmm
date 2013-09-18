classdef OctreeNode < handle
    
    properties
    
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
        
        % this node index
        Index
        
        % this node's lists 
        InteractionList
        NearFieldList
        
        
        %%%%%%%%%%%% STATS %%%%%%%%%%%%%%%%%%%%%%
        
        % the outgoing representation
        PsiVector
        
        % the incoming representation
        PhiVector
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        InterpolationPoints
        
        ProjMatrix
        EvalMatrix

        OutgoingSkeleton
        OutgoingSkeletonSize
        
        IncomingSkeleton
        IncomingSkeletonSize
        
        OIMatrices
        
        % for now, assuming that in a leaf the outgoing and incoming
        % skeletons are just all of the points in the leaf
    
    end
    
    methods
        
        function obj = OctreeNode(begin, count, depth, min_val, max_val, index)
            
            if (nargin > 0)
                if (count > 0) 
                    obj.Begin = begin;
                    obj.Count = count;
                    obj.Depth = depth;
                    obj.MinVal = min_val;
                    obj.MaxVal = max_val;
                    obj.End = begin + count - 1;
                    obj.HasLeft = false;
                    obj.HasRight = false;
                    obj.Index = index;
                    
                    %obj.InteractionList = OctreeNode.empty(2, 0);
                    obj.InteractionList = int32.empty(2,0);
                    
                    %obj.NearFieldList = OctreeNode.empty(2, 0);
                    obj.NearFieldList = int32.empty(2,0);
                    
                else
                    obj.Index = -1;
                end
            else
                obj.Index = -1;
            end
        end
                
        function CopySkeletons(this, child)
           
            this.ProjMatrix = child.ProjMatrix;
            this.OutgoingSkeleton = child.OutgoingSkeleton;
            this.OutgoingSkeletonSize = child.OutoingSkeletonSize;
            this.IncomingSkeleton = child.IncomingSkeleton;
            this.IncomingSkeletonSize = child.IncomingSkeletonSize;
           
            this.InterpolationPoints = child.InterpolationPoints;
            
        end
        
        function res = is_empty(this)
           
            if (this.Index > 0)
                res = false;
            else
                res = true;
            end
            
        end
        
        function Print(this)
           
            'Node:'
            begin = this.Begin
            count = this.Count
            
        end
        
        function res = IsLeaf(this)
           
            if (this.is_empty())
                res = false;
            else 
                res = ~(this.HasLeft || this.HasRight);
            end
            
        end
        
        function res = WellSeparated(this, other)
           
            % I'm assuming that if they don't touch, then they're
            % well-separated
            % I think this assumes that they're on the same level of the
            % tree
            
            res = false;
            eps = 1e-16;
            
            % make sure they're not the same node and that they are on the
            % same level of the tree
            if (this.Index ~= other.Index && this.Depth == other.Depth)
                
                if (abs(this.MaxVal - other.MinVal) > eps && abs(other.MaxVal - this.MinVal) > eps) 
                    res = true;
                end
                
            end
            
        end
        
        function dist = MinDistance(this, point)
            
            lower = this.MinVal - point;
            upper = point - this.MaxVal;
            
            dist = 0.5 * (lower + abs(lower) + upper + abs(upper));
            
        end
        
                
    end
    
end

