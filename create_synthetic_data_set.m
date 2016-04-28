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

% Henry Schuetze 
% Institute for Neuro- and Bioinformatics
% University of Luebeck, Germany
% Henry.Schuetze@uni-luebeck.de
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