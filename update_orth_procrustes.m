% params = update_orth_procrustes(params)
%
% Performs an orthogonal dictionary update by solving the orthogonal 
% procrustes problem subject to data matrix and sparse coefficient matrix.
%
% INPUT: params
% ======
% 
% params.X (required):
%   data matrix (num_dims x num_samples)
%
% params.U (required):
%   orthogonal dictionary (num_dims x num_atoms)
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
function params = update_orth_procrustes(params)

params.A = sparse_coefficients(params.X, params.U, ...
    params.sparse_mode, params.sparsity_param);
params.U = solve_orth_procrustes(params.X * params.A');