function [ error, ran ] = BoxCompression(A, epsilon, filter)

    disp('Decomposing matrix');
    % this is for the uniform sampler
    if (numel(filter) == numel(A))
        Asub = A .* filter;
    else
        % assumming the row sampler
        % we assume that the rows in filter are sorted here, I think, but
        % it shouldn't matter
        Asub = A(filter,:);
        
    end
    
    [proj, skeleton] = InterpolativeDecomposition(Asub, epsilon);
    
    ran = numel(skeleton);
    
    disp('Reconstructing approximation');
    Aapprox = A(:, skeleton) * proj;

    error = norm(A - Aapprox) / norm(A);
    
end

