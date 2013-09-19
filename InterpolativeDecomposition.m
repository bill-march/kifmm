function [Acol, proj, skeleton] = InterpolativeDecomposition(A)
%InterpolativeDecompostion -- Decompose the matrix as a linear combination
%of k of its columns
% A col is a matrix consisting of k columns of A
% proj is a k x N matrix (containing the k x k identity matrix as a
% submatrix
% skeleton is an array of the indices of the columns of A that appear in
% Acol

% IMPORTANT: we're assuming that the columns that appear in Acol are in the
% same order they were in in A

    %k = rank(A, 1e-15);
    

    % As a debugging measure, compute the trivial decomposition that
    % consists of all columns
    Acol = A;
    proj = eye(size(A));
    skeleton = 1:size(A,2);
    
    

end

