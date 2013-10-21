function [ errors, ranks ] = CompressionTest2( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    box_width = (leaf_size * gamma(dimensions/2 + 1) / (num_data * pi^(dimensions/2)))^(1/dimensions);


    % Just uniform donut inside the far field
    
    data = rand(dimensions,num_data);
    %data = sort(data);

    % the size of a ball that contains leaf_size points in expectation
    
    leaf_corner = max(rand(dimensions, 1) - box_width,0);
    box_center = leaf_corner + box_width / 2;
    
    
    [my_box, far_field] = FindPointsInBox(data, box_center, box_width);
    
    %my_box = 1:leaf_size;
    %far_field = 3*leaf_size:num_data;
    
    % Grow the box until it contains enough elements
    while (numel(my_box) < leaf_size)
        
        box_step = 1e-5;
        
        box_width = box_width + box_step;
        leaf_corner = max(rand(dimensions, 1) - box_width,0);
        box_center = leaf_corner + box_width / 2;

        [my_box, far_field] = FindPointsInBox(data, box_center, box_width);        
        
    end

    
    
    
    %Better separated points
    

    %uniformly distributed near field, scaled to live in small box
    
    
    near_data = rand(dimensions, leaf_size);
    near_data = near_data * box_width;
    
    far_data = rand(dimensions, num_data);
    far_trans = zeros(dimensions,1);
    far_trans(1) = 10 * box_width;
    far_data = far_data * 2 * box_width;
    far_data = bsxfun(@plus, far_data, far_trans);
    
    
    data = [near_data, far_data];
    my_box = 1:leaf_size;
    far_field = leaf_size+1:leaf_size+num_data;

    
    
    errors = zeros(numel(epsilons), numel(sampling_probs));
    ranks = zeros(numel(epsilons), numel(sampling_probs));
    
    disp('Building full matrix');
    A = kernel.eval_mat(data, far_field, my_box);
    
    for i = 1:numel(epsilons)
       
        for j = 1:numel(sampling_probs)

            sampler = RowSampler(sampling_probs(j));
            
            disp(['Doing ', num2str(i), ', ', num2str(j)]);
            [errors(i,j), ranks(i,j)] = BoxCompression(A, epsilons(i), sampler.GenFilter(1:numel(far_field)));

        end
        
    end
    
end




end

