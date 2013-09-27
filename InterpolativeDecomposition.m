function [proj, skeleton] = InterpolativeDecomposition(A, k)
%InterpolativeDecompostion -- Decompose the matrix as a linear combination
%of k of its columns
% proj is a k x N matrix (containing the k x k identity matrix as a
% submatrix
% skeleton is an array of the indices of the columns of A that appear in
% Acol
% A is the matrix to be decomposed
% k is the number of columns to use in the decomposition

    % for now, just going with the bare QR decomposition version
    
    % we can't take more columns of A than exist
    r = rank(A);
    k = min(r,k);

    [~, R, perm] = qr(A,0);
    
    % I think R should be in decreasing order of diagonal values
    % Now, R should be some linear combination of its own columns
    % These same columns will form the ID of A
    
    skeleton = perm(1:k);
    
    % Want X s.t. R(:,1:k) X = R
    size(A)
    proj = R(:,1:k) \ R
    
    proj(:,perm) = proj;
    
    %Aapprox = A(:,skeleton) * proj;
    
    %err = norm(A - Aapprox)
    
    %if (err > 1e-15)
    %    'found it'
    %end
    
end

