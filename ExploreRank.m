function [ S ] = ExploreRank( num_data, dimensions, cube_size, leaf_size, kernel )

    data = rand(dimensions, num_data) .* cube_size;

    box_center = ones(dimensions,1) .* cube_size / 2;
    box_width = (leaf_size * gamma(dimensions/2 + 1) / num_data)^(1/dimensions) * pi^(-1/2) * cube_size;
    
    
    [near_points, far_points] = FindPointsInBox(data, box_center, box_width);
    
    K = kernel.eval_mat(data, far_points, near_points);
    
    S = svd(K);

end

