% params = update_GFOSC_line_search(params)
%
% Performs an orthogonal dictionary update according to the
% Geodesic-Flow OSC algorithm using a backtracking line search.
%
% [1] Schütze, H., Barth, E., and Martinetz, T. "Learning Orthogonal Sparse
% Representations by Using Geodesic Flow Optimization: Neural Networks 
% (IJCNN), 2015 International Joint Conference on, 15540:1-8
%
% [2] https://en.wikipedia.org/wiki/Backtracking_line_search
%
% INPUT: params
% ======
% 
% params.x (required):
%   data sample (num_dims x 1)
%
% params.U (required):
%   complete orthogonal dictionary (num_dims x num_dims)
%
% params.alpha0 (required, params.alpha0 > 0):
%   maximum candidate step size value of backtracking line search
%
% params.tau (required, 0 < params.tau < 1):
%   search control parameter of backtracking line search
%
% params.c (required, , 0 < params.c < 1):
%   search control parameter of backtracking line search
%
% params.sparse_mode (required):
%   a string with value either 'column_k-sparse' or 'hard_thresh' selecting
%   the sparse model
%
% params.sparsity_param (required):
%   sparsity parameter corresponding to the sparse model, i.e., either the 
%   number of non-zero coefficients (if sparse_mode == 'column_k-sparse') 
%   or the hard threshold (if sparse_mode == 'hard_thresh')
%
% OUTPUT: params
% =======
%
% params.U:
%   updated dictionary

% Henry Schuetze 
% Institute for Neuro- and Bioinformatics
% University of Luebeck, Germany
% Henry.Schuetze@uni-luebeck.de
function params = update_GFOSC_line_search(params)

alpha = params.alpha0;

a_hat = sparse_coefficients(params.x, params.U, params.sparse_mode, params.sparsity_param);

x_hat = params.U*a_hat;
grad_E_A = x_hat*params.x' - params.x*x_hat';

t = -params.c*(-grad_E_A(:)' * grad_E_A(:));
cost_at_U = cost_function(params.x, params.U, params.sparse_mode, params.sparsity_param);

delta_U = expm(alpha*(-grad_E_A));
cost_at_delta_UU = cost_function(params.x, delta_U*params.U, params.sparse_mode, params.sparsity_param);

while (cost_at_U - cost_at_delta_UU) < alpha*t
    alpha = params.tau*alpha;
    delta_U = expm(-alpha*grad_E_A);
    cost_at_delta_UU = cost_function(params.x, delta_U*params.U, ...
        params.sparse_mode, params.sparsity_param);
end

params.U = delta_U*params.U;
