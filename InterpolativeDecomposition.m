function [Acol, proj, skeleton] = InterpolativeDecomposition(A)
%InterpolativeDecompostion -- Decompose the matrix as a linear combination
%of k of its columns
% A col is a matrix consisting of k columns of A
% proj is a k x N matrix (containing the k x k identity matrix as a
% submatrix
% skeleton is an array of the indices of the columns of A that appear in
% Acol

    k = rank(A, 1e-15);
    
    % now, we need ell Gaussian iid columns, but how to choose ell?

    
    
    

end

