% params = update_lambdaOSC_fast(params)
%
% Performs an orthogonal dictionary update similiar to the OSC algorithm.
% This function implements a naive straight-forward extension of OSC to
% meet the unconstrained Lagrangian model.
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
% params.sparsity_param (required):
%   the hard threshold (independent of params.sparse_mode)
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