% params = update_GFOSC_cooling_learnrate(params)
%
% Performs an orthogonal dictionary update according to the
% Geodesic-Flow OSC algorithm using a cooling learn rate.
%
% [1] Schütze, H., Barth, E., and Martinetz, T. "Learning Orthogonal Sparse
% Representations by Using Geodesic Flow Optimization: Neural Networks 
% (IJCNN), 2015 International Joint Conference on, 15540:1-8
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
function params = update_GFOSC_cooling_learnrate(params)

eps_t = params.eps_i * (params.eps_f/params.eps_i)^(params.t/params.t_max);

a_hat = sparse_coefficients(params.x, params.U, params.sparse_mode, ...
    params.sparsity_param);
x_hat = params.U*a_hat;

grad_E_A = x_hat*params.x' - params.x*x_hat';
delta_U = expm(-eps_t*grad_E_A);

params.U = delta_U*params.U;