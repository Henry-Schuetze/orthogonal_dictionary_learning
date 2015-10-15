% [X, A] = create_synthetic_data_set(C, num_samples, k, synth_mode)
%
% Creates a synthetic dataset which has a sparse representation with
% respect to the dictionary C
%
% INPUT:
% ======
% C (required):
%   dictionary (num_dims x num_atoms)
%
% num_samples (required):
%   the number of samples (columns of data matrix X) to be generated
%
% k (required):
%   sparsity parameter (if synth_mode == 'column_k-sparse' it defines the
%   exact number of non-zero coefficients, if synth_mode ==
%   'bernoulli_gaussian' it defines the probability of a coefficient beeing
%   non-zero as prob_bernoulli = k/num_atoms)
%
% synth_mode (required):
%   a string with value either 'column_k-sparse' or 'bernoulli_gaussian' 
%   selecting the synthesis model
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
function [X, A] = create_synthetic_data_set(C, num_samples, k, synth_mode)

num_atoms = size(C, 2);

Z = zeros(num_atoms, num_samples);
switch synth_mode
    case 'column_k-sparse'
        for i = 1:num_samples
            idx_vec = randperm(num_atoms, k);
            Z(idx_vec, i) = 1;
        end
        
    case 'bernoulli_gaussian'
        prob_bernoulli = k/num_atoms;
        Z(rand(num_atoms, num_samples) <= prob_bernoulli) = 1;
end

A = randn(num_atoms, num_samples);
A = A.*Z;

X = C*A;