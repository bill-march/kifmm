function [proj, skeleton] = InterpolativeDecomposition(A, k)
%InterpolativeDecompostion -- Decompose the matrix as a linear combination
%of k of its columns
% proj is a k x N matrix (containing the k x k identity matrix as a
% submatrix
% skeleton is an array of the indices of the columns of A that appear in
% Acol
% A is the matrix to be decomposed
% k is the number of columns to use in the decomposition


    % this is the trivial version, in which we take all of the columns of A
    proj = eye(size(A,2));
    skeleton = 1:size(A,2);
    
    % for now, just going with the bare QR decomposition version
    
    [Q, R, perm] = qr(A,0);
    
    
    
    
    

end

