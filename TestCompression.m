function [ errors, ranks ] = TestCompression( num_data, dimensions, leaf_size, kernel, epsilons, sampling_probs )

    rng('shuffle');

%     data_scale = 1;
% 
% 
%     % Just uniform donut inside the far field
%     
%     data = rand(dimensions,num_data) .* data_scale;
%     %data = sort(data);
%  
%     % the size of a ball that contains leaf_size points in expectation
%     
%     %leaf_corner = max(data_scale .* rand(dimensions, 1) - box_width,0);
%     %box_center = leaf_corner + box_width / 2;
%     
%     % just put the leaf in the center
%     box_center = ones(dimensions,1) * data_scale / 2;
%     
%     
%     [my_box, far_field] = FindPointsInBox(data, box_center, box_width);
%     
%     %my_box = 1:leaf_size;
%     %far_field = 3*leaf_size:num_data;
%     
%     % Grow the box until it contains enough elements
%     while (numel(my_box) < leaf_size)
%         
%         'expanding box'
%         
%         box_step = 1e-5;
%         
%         box_width = box_width + box_step;
%         
%         [my_box, far_field] = FindPointsInBox(data, box_center, box_width);        
%         
%     end
% 
%     disp(box_width);
%     
%     
%     %Better separated points
%     

    %uniformly distributed near field, scaled to live in small box
    
    
%     near_data = rand(dimensions, leaf_size);
%     near_data = near_data * box_width;
%     
%     far_data = rand(dimensions, num_data);
%     far_trans = zeros(dimensions,1);
%     far_trans(1) = 10 * box_width;
%     far_data = far_data * 2 * box_width;
%     far_data = bsxfun(@plus, far_data, far_trans);
%     
%     
%     data = [near_data, far_data];
%     my_box = 1:leaf_size;
%     far_field = leaf_size+1:leaf_size+num_data;
% 
%     
    
    ws_index = 2;


    % Generating points in the sphere
    
    %sphere_width = (leaf_size * gamma(dimensions/2 + 1) / (num_data * pi^(dimensions/2)))^(1/dimensions) * data_scale;
    sphere_width = 1;
    
    near_points = randsphere(dimensions, leaf_size, sphere_width);

    data_scale = (num_data / leaf_size)^(1/dimensions) * sphere_width
    
    far_points = randsphere(dimensions, num_data, data_scale);
    data = [near_points, far_points];
    
    [~, far_field] = FindPointsInBox(data, zeros(dimensions,1), sphere_width, ws_index);
    my_box = 1:leaf_size;
    disp(['Near points: ', num2str(numel(my_box))]);
    disp(['Far points: ', num2str(numel(far_field))]);
    
    errors = zeros(numel(epsilons), numel(sampling_probs));
    ranks = zeros(numel(epsilons), numel(sampling_probs));
    
    disp('Building full matrix');
    A = kernel.eval_mat(data, far_field, my_box);

    max(max(A))
    
    for i = 1:numel(epsilons)
       
        for j = 1:numel(sampling_probs)

            sampler = RowSampler(sampling_probs(j));
            
            disp(['Doing ', num2str(i), ', ', num2str(j)]);
            [errors(i,j), ranks(i,j)] = BoxCompression(A, epsilons(i), sampler.GenFilter(1:numel(far_field)));

        end
        
    end
    
end



