function [proj, skeleton] = InterpolativeDecomposition(A, eps)
%InterpolativeDecompostion -- Decompose the matrix as a linear combination
%of k of its columns
% proj is a k x N matrix (containing the k x k identity matrix as a
% submatrix
% skeleton is an array of the indices of the columns of A that appear in
% Acol
% A is the matrix to be decomposed
% k is the number of columns to use in the decomposition

    % for now, just going with the bare QR decomposition version
    k = rank(A,eps);

    [~, R, perm] = qr(A,0);
    
    % I think R should be in decreasing order of diagonal values
    % Now, R should be some linear combination of its own columns
    % These same columns will form the ID of A
    
    skeleton = perm(1:k);
    
    B(:,perm) = R;
    
    % Want X s.t. B(:,1:k) X = B
    proj = B(:,skeleton) \ B;
    
    Aapprox = A(:,skeleton) * proj;
            
end

