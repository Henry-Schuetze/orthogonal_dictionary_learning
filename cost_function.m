% cost_val = cost_function(X, U, sparse_mode, sparsity_param)
%
% Evaluates the cost function for orthogonal sparse coding. The sparse
% coefficient matrix A is computed by the function sparse_coefficients(X,
% U, sparse_mode, sparsity_param). If sparse_mode == 'column_k-sparse' then
% the cost function ||X - U*A||²_F is evaluated. If sparse_mode ==
% 'hard_thresh' then the cost function ||X - U*A||²_F + lambda*||A||_0 is
% evaluated (where each entry of A with absolute value less than 1e-4 is 
% considered to be zero).

%
% INPUT:
% ======
% X (required):
%   data matrix (num_dims x num_samples)
%
% U (required):
%   orthogonal dictionary U (num_dims x num_atoms)
%
% sparse_mode (required):
%   a string with value either 'column_k-sparse' or 'hard_thresh' selecting
%   the sparse model
%
% sparsity_param (required):
%   sparsity parameter corresponding to the sparse model, i.e., either k,
%   the number of non-zero coefficients (if sparse_mode == 
%   'column_k-sparse') or lambda, the hard threshold (if sparse_mode == 
%   'hard_thresh').
%
% OUTPUT:
% =======
% cost_val
%   the cost function value

% Henry Schuetze 
% Institute for Neuro- and Bioinformatics
% University of Luebeck, Germany
% Henry.Schuetze@uni-luebeck.de
function cost_val = cost_function(X, U, sparse_mode, sparsity_param)

% compute the sparse coefficient matrix A for data matrix X with respect to
% dictionary U
A = sparse_coefficients(X, U, sparse_mode, sparsity_param);

% compute the approximation error term ||X - U*A||²_F
cost_val = sum(sum((X - U*A).^2));

if strcmp(sparse_mode, 'hard_thresh')
    % add lambda*||A||_0 to the approximation error term
    cost_val = cost_val + sparsity_param*sum(sum(abs(A) >= 1e-14));
end