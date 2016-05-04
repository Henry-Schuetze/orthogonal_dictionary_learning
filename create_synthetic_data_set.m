% [X, A] = create_synthetic_data_set(C, num_samples, sparsity_param, ...
%   synthetic_mode)
%
% Creates a synthetic dataset of sparse samples, where each sample is 
% linearly combined by some small subset of atoms of the (not necessarily 
% orthogonal) dictionary C. The locations of non-zero coefficients are
% selected uniformly at random. The non-zero coefficients are drawn (i.i.d)
% from a standard Gaussian distribution.
%
% INPUT:
% ======
% C (required):
%   dictionary (num_dims x num_atoms)
%
% num_samples (required):
%   number of desired samples (columns of data matrix X) to be generated
%
% sparsity_param (required):
%   sparsity parameter corresponding to the synthetic mode (if
%   synthetic_mode == 'column_k-sparse' it defines for each sample the 
%   exact number of non-zero coefficients, if synthetic_mode ==
%   'bernoulli_gaussian' it defines the probability of a coefficient to be
%   non-zero which is sparsity_param/num_atoms)
%
% synthetic_mode (required):
%   a string with value either 'column_k-sparse' or 'bernoulli_gaussian' 
%   selecting the synthetic mode
%
% OUTPUT:
% =======
% X:
%   data matrix (num_dims x num_samples)
%
% A:
%   sparse coefficient matrix (num_atoms x num_samples) by which the data 
%   matrix X is synthesized

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
function [X, A] = create_synthetic_data_set(C, num_samples, ...
    sparsity_param, synthetic_mode)

num_atoms = size(C, 2);

Z = zeros(num_atoms, num_samples);
switch synthetic_mode
    case 'column_k-sparse'
        for i = 1:num_samples
            idx_vec = randperm(num_atoms, sparsity_param);
            Z(idx_vec, i) = 1;
        end
        
    case 'bernoulli_gaussian'
        prob_bernoulli = sparsity_param/num_atoms;
        Z(rand(num_atoms, num_samples) <= prob_bernoulli) = 1;
end

A = randn(num_atoms, num_samples);
A = A.*Z;

X = C*A;