% Main function for kernel-independent fmm
function [] = kifmm(Sources, Targets, Kernel)



    Souce_Tree = BuildTree(Sources);
    Target_Tree = BuildTree(Targets);
    
    

    
end

function points = SampleFarField()

    

end

function PrecomputeRepresentation(tree, node, kernel)
           % only gets called for leaves
           
           % First, sample some points from this node's far-field
           interpolation_points = SampleFarField(tree, node, num_interpolation_points);
            
           direct_matrix = zeros(size(interpolation_points, 2), node.Count);
           direct_matrix_trans = zeros(node.Count, size(interpolation_points, 2));
           
           for i = 1:num_interpolation_points
              
               for j = node.Begin:node.End
                  
                   direct_matrix(i,j) = kernel.eval(data(interpolation_points(i)), data(j));
                   direct_matrix_trans(j,i) = kernel.eval(data(j), data(interpolation_points(i)));
                   
               end
               
           end
           
           [Acol proj] = InterpolativeDecomposition(direct_matrix);
           
           % For now, just letting the outgoing representation be the
           % points sampled above
           % So, the outgoing skeleton is all points in the node, since
           % we're assuming that the leaf nodes contain few points
           
           node.ProjMatrix = proj;
           node.Brow = direct_matrix_trans;
           
           % applying the above assumption to the incoming representation,
           % we have that 
           
           
end




