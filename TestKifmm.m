function [rms_error, max_error, times] = TestKifmm(num_data, leaf_size, tree_depth, skel_sizes)

% outputs average and max squared error (relative to naive algorithm) for
% each value in skel_sizes and for given fixed tree parameters

    % uniformly random data in the unit interval
    data = rand(1,num_data);
    % unit charges for now
    charges = ones(1,num_data);
    
    kernel_band = 1.0;    
    kernel = GaussianKernel(kernel_band);
    
    
    num_skels = size(skel_sizes,2);
    
    rms_error = zeros(1,num_skels);
    max_error = zeros(1,num_skels);
    max_rel_error = zeros(1,num_skels);
    times = zeros(2,num_skels);
    
    % minimum size of a node in the tree
    min_node_size = 1e-8;
    
    for i = 1:num_skels
        
        fprintf('skeleton %d\n', i)

        tic;
        kifmm = KIFMM(data, leaf_size, min_node_size, tree_depth, skel_sizes(i), kernel);
        tree_pot = kifmm.ComputePotentials(charges);
        times(1,i) = toc;
        
        tic;
        naive_pot = kifmm.ComputePotentialsNaive(data, charges);
        times(2,i) = toc
        
        rms_error(i) = sqrt((tree_pot - naive_pot)' * (tree_pot - naive_pot) / (naive_pot' * naive_pot));
        max_error(i) = max(abs(tree_pot - naive_pot) * num_data / sum(abs(naive_pot)));
        
    end



end

