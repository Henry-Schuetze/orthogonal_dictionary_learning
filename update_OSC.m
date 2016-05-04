% params = update_OSC(params)
%
% Performs an orthogonal dictionary update according to the OSC algorithm.
% This function also implements a naive straight-forward extension to meet
% the unconstrained Lagrangian model.
%
% Consider to use update_OSC_fast.mexa64 or update_lambdaOSC_fast.mexa64 
% instead of update_OSC.m for improved runtime performance.
%
% [1] Schütze, H., Barth, E., and Martinetz, T., "Learning Efficient Data
% Representations with Orthogonal Sparse Coding", IEEE Transactions on 
% Computational Imaging, 2016 (accepted)
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

% Copyright © 2016 Henry Schuetze
% Institute for Neuro- and Bioinformatics
% University of Luebeck, Germany
% Henry.Schuetze@uni-luebeck.de
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
function params = update_OSC(params)

[a, seq_vec] = sparse_coefficients(params.x, params.U, ...
    params.sparse_mode, params.sparsity_param);
eps_t = params.eps_i * (params.eps_f/params.eps_i)^(params.t/params.t_max);

for k = 1 : size(params.U,2)
    for l = 1 : (k-1)
        % orthogonalize atom u_k subject to updated u_1,...,u_(k-1)
        params.U(:, seq_vec(k)) = params.U(:, seq_vec(k)) - ...
            (params.U(:, seq_vec(k))'*params.U(:, seq_vec(l))) * params.U(:,seq_vec(l));
    end
    
    if (strcmp(params.sparse_mode, 'column_k-sparse') && k <= params.sparsity_param) || ...
            (strcmp(params.sparse_mode, 'hard_thresh') && abs(a(seq_vec(k))) > params.sparsity_param)
        % Hebbian-like main update
        params.U(:, seq_vec(k)) = params.U(:, seq_vec(k)) + ...
            eps_t * (params.U(:, seq_vec(k))'* params.x) * params.x;
    end
    
    % normalize atom u_k
    params.U(:, seq_vec(k)) = params.U(:, seq_vec(k)) / sqrt(sum(params.U(:, seq_vec(k)).^2));
    
    % update residual
    params.x = params.x - (params.U(:, seq_vec(k))'*params.x)*params.U(:, seq_vec(k));
end