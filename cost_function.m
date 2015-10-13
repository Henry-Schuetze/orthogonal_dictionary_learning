% cost_val = cost_function(params)
%
% Evaluates the cost function for orthogonal sparse coding
%
% INPUT: params
% ======
% X (required):
%   data matrix (num_dims x num_samples)
%
% U (required):
%   orthogonal dictionary U (num_dims x num_atoms)
%
% sparse_mode (required):
%   either 'column_k-sparse' or 'hard_thresh'
%
% params.k or params.thresh (one of both required):
%   sparsity parameter depending on params.sparse_mode
%
% OUTPUT:
% =======
% cost_val
%   cost function value
function cost_val = cost_function(X, U, sparse_mode, sparsity_param)

% compute sparse coefficient matrix A of X w.r.t. U
A = sparse_coefficients(X, U, sparse_mode, sparsity_param);

% compute ||X - U*A||²_F
cost_val = sum(sum((X - U*A).^2));

if strcmp(sparse_mode, 'hard_thresh')
    % add lambda*||A||_0 to cost_val
    cost_val = cost_val + sparsity_param*sum(sum(abs(A) >= 1e-14));
end