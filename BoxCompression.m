function [ error, ran ] = BoxCompression(A, epsilon, filter)

    disp('Decomposing matrix');
    Asub = A .* filter;
    [proj, skeleton] = InterpolativeDecomposition(Asub, epsilon);
    
    ran = numel(skeleton);
    
    disp('Reconstructing approximation');
    Aapprox = Asub(:, skeleton) * proj;

    error = norm(A - Aapprox);
    
end

