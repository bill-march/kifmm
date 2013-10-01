function [rms_error, max_error, times, num_evals, naive_time, naive_num_evals] = TestKifmm(num_data, leaf_size, tree_depth, epsilons, probs)

% outputs average and max squared error (relative to naive algorithm) for
% each value in skel_sizes and for given fixed tree parameters

    % uniformly random data in the unit interval
    data = rand(1,num_data);
    % unit charges for now
    charges = ones(1,num_data);
    
    kernel_band = 1.0;    
    kernel = GaussianKernel(kernel_band);
    
    
    num_epsilons = size(epsilons,2);
    num_probs = size(probs,2);
    
    rms_error = zeros(num_epsilons, num_probs);
    max_error = zeros(num_epsilons, num_probs);
    times = zeros(num_epsilons, num_probs);
    num_evals = zeros(num_epsilons, num_probs);
    
    % minimum size of a node in the tree
    min_node_size = 1e-8;
    
    
    sampler = UniformSampler(0.1);
    kifmm = KIFMM(data, leaf_size, min_node_size, tree_depth, 1e-5, kernel, sampler);
    tic;
    naive_pot = kifmm.ComputePotentialsNaive(data, charges);
    naive_time = toc;
    
    % we don't assume anything about symmetry
    naive_num_evals = num_data * num_data;
            
    for i = 1:num_epsilons
        
        for j = 1:num_probs
        
            fprintf('epsilon %g, s %g\n', epsilons(i), probs(j))

            tic;
            sampler = UniformSampler(probs(j));
            kifmm = KIFMM(data, leaf_size, min_node_size, tree_depth, epsilons(i), kernel, sampler);
            tree_pot = kifmm.ComputePotentials(charges);
            times(i,j) = toc;

            rms_error(i,j) = sqrt((tree_pot - naive_pot)' * (tree_pot - naive_pot) / (naive_pot' * naive_pot));
            max_error(i,j) = max(abs(tree_pot - naive_pot) * num_data / sum(abs(naive_pot)));
            num_evals(i,j) = kifmm.NumKernelEvaluations;
            
        end
        
    end

end

