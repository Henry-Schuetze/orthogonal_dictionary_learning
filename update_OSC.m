% params = update_OSC(params)
%
% Performs an orthogonal dictionary update according to the OSC algorithm.
%
% [1] Schütze, H., Barth, E., and Martinetz, T. (2016), "Learning Efficient 
% Data Representations with Orthogonal Sparse Coding", IEEE Transactions on 
% Computational Imaging
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
% params.eps_i (required):
%   initial learning rate
%
% params.eps_f (required):
%   final learning rate
%
% params.t (required):
%   iteration counter of performed learning steps
%
% params.params.t_max (required):
%   maximal number of learning steps
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
function params = update_OSC(params)

[~, seq_vec] = sparse_coefficients(params.x, params.U, ...
    params.sparse_mode, params.sparsity_param);
eps_t = params.eps_i * (params.eps_f/params.eps_i)^(params.t/params.t_max);

for k = 1 : size(params.U,2)
    for l = 1 : (k-1)
        % orthogonalize atom u_k subject to update u_1,...,u_(k-1)
        params.U(:, seq_vec(k)) = params.U(:, seq_vec(k)) - ...
            (params.U(:, seq_vec(k))'*params.U(:, seq_vec(l))) * params.U(:,seq_vec(l));
    end
    
    if k <= params.sparsity_param
        % Hebbian-like main update
        params.U(:, seq_vec(k)) = params.U(:, seq_vec(k)) + eps_t * (params.U(:, seq_vec(k))'* params.x) * params.x;
    end
    
    % normalize atom u_k
    params.U(:, seq_vec(k)) = params.U(:, seq_vec(k)) / sqrt(sum(params.U(:, seq_vec(k)).^2));
    
    % update residual
    params.x = params.x - (params.U(:, seq_vec(k))'*params.x)*params.U(:, seq_vec(k));
end