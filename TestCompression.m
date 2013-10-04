function [ errors, ranks ] = TestCompression( num_data, dimensions, leaf_size, kernel, epsilons, sampling_probs )

    data = rand(dimensions,num_data);
    %data = sort(data);

    % the size of a ball that contains leaf_size points in expectation
    box_width = (leaf_size * gamma(dimensions/2 + 1) / (num_data * pi^(dimensions/2)))^(1/dimensions);
    
    leaf_corner = max(rand(dimensions, 1) - box_width,0);
    box_center = leaf_corner + box_width / 2;

    
    
    [my_box, far_field] = FindPointsInBox(data, box_center, box_width);
    
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

