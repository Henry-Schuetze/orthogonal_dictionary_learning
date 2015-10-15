% [A, I_sorted] = sparse_coefficients(X, U, sparse_mode, sparsity_param)
%
% Computes the sparse coefficients for data matrix X (samples are stored
% columnwise) subject to the orthogonal dictionary U.
%
% INPUT:
% ======
% X (required):
%   data matrix (num_dims x num_samples) or single sample (num_dims x 1)
%
% U (required):
%   orthogonal dictionary U (num_dims x num_atoms)
%
% sparse_mode (required):
%   a string with value either 'column_k-sparse', 'hard_thresh' or 
%   'fix_support'selecting the sparse model
%
% sparsity_param (required):
%   sparsity parameter corresponding to the sparse model, i.e., either the 
%   number of non-zero coefficients (if sparse_mode == 'column_k-sparse') 
%   or the hard threshold (if sparse_mode == 'hard_thresh')
%
% OUTPUT:
% =======
% A:
%   sparse coefficients
%
% I_sorted:
%   indices which sort the "plain" coefficients (i.e. U'X) columnwise in 
%   descending order (caution: only set when sparse_mode == 
%   'column_k-sparse')

% Henry Schuetze 
% Institute for Neuro- and Bioinformatics
% University of Luebeck, Germany
% Henry.Schuetze@uni-luebeck.de
function [A, I_sorted] = sparse_coefficients(X, U, sparse_mode, sparsity_param)

num_samples = size(X, 2);
dict_size = size(U, 2);

A = U'*X;

switch sparse_mode
    case 'column_k-sparse'
        [~, I_sorted] = sort(abs(A), 'descend');
        
        if num_samples > 1
            % X is a matrix
            %
            % fast implementation to set the (num_dims-k) smallest
            % |coefficients| to zero for each column in A:
            row_coord = uint32(reshape(I_sorted((sparsity_param+1):dict_size,:), ...
                (dict_size-sparsity_param)*num_samples, 1));
            col_coord = uint32(kron((1:num_samples)', ...
                ones(dict_size-sparsity_param,1)));
            I_zero = uint32(sub2ind([dict_size, num_samples], ...
                row_coord, col_coord));
            A(I_zero) = 0;

        else
            % X is a vector (single sample)
            A(I_sorted((sparsity_param+1):end)) = 0;
        end
        
    case 'fix_support'
        I_zero = params.Z ~= 0;
        A(I_zero) = 0;
    
    case 'hard_thresh',
        I_zero = abs(A) < sparsity_param;
        A(I_zero) = 0;
end