function [ errors, ranks ] = TestCompression( num_data, leaf_size, kernel, epsilons, sampling_probs )

    data = rand(1,num_data);
    data = sort(data);
    
    start_ind = randi(num_data-leaf_size,1);
    
    my_box = start_ind:start_ind+leaf_size-1;
    my_box_size = data(start_ind+leaf_size-1) - data(start_ind);
    
    % far-field is more than a box width away
    far_field_l = find(data < data(start_ind) - my_box_size);
    far_field_r = find(data > data(start_ind) + 2 * my_box_size);

    far_field = [far_field_l, far_field_r];
    
    errors = zeros(numel(epsilons), numel(sampling_probs));
    ranks = zeros(numel(epsilons), numel(sampling_probs));
    
    disp('Building full matrix');
    A = kernel.eval_mat(data, far_field, my_box);

    for i = 1:numel(epsilons)
       
        for j = 1:numel(sampling_probs)

            disp(['Doing ', num2str(i), ', ', num2str(j)]);
            [errors(i,j), ranks(i,j)] = BoxCompression(A, epsilons(i), sampling_probs(j));

        end
        
    end
    
    
end

