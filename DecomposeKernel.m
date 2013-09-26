function [Acol, proj, skeleton] = DecomposeKernel(data, row_inds, col_inds, kernel)
%DecomposeKernel -- Decompose the matrix as a linear combination
%of k of its columns
% A col is a matrix consisting of k columns of A
% proj is a k x N matrix (containing the k x k identity matrix as a
% submatrix
% skeleton is an array of the indices of the columns of A that appear in
% Acol
% data is the data  matrix
% row and col_inds contain the indices we are decomposing the submatrix of
% kernel is a class with an eval() method that takes two data points and
% outputs the corresponding entry of the Gram matrix

% IMPORTANT: we're assuming that the columns that appear in Acol are in the
% same order they were in in A

    %k = rank(A, 1e-15);
    

    % As a debugging measure, compute the trivial decomposition that
    % consists of all columns
    % so, k is the number of columns of A
    proj = eye(size(col_inds,2));
    skeleton = col_inds;
    
    Acol = zeros(size(row_inds,2), size(col_inds, 2));
    
    for i = 1:size(row_inds, 2)
       
        for j = 1:size(col_inds, 2)
           
            Acol(i,j) = kernel.eval(data(row_inds(i)), data(col_inds(j)));
            
        end
        
    end

end

