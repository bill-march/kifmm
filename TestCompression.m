function [ errors, ranks ] = TestCompression( num_data, dimensions, leaf_size, kernel, epsilons, sampling_probs )

    data = rand(dimensions,num_data);
    %data = sort(data);
    
    leaf_width = (leaf_size / num_data)^(1/dimensions) * gamma(dimensions/2) / pi^(dimensions/2)
    leaf_corner = max(rand(dimensions, 1) - leaf_width,0);
    leaf_center = leaf_corner + leaf_width / 2;
    
    % now, find the points that live in the box
    [my_box, far_field] = FindPointsInBox(data, leaf_center, leaf_width);
    
    disp(numel(my_box))
    
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

