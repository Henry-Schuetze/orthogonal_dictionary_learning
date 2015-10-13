function [X, Z] = create_synthetic_data_set(U, num_samples, k, synth_mode)

num_atoms = size(U, 2);
assert(k <= num_atoms);

switch synth_mode
    case 'column_k-sparse'
        Z = zeros(num_atoms, num_samples);

        for i = 1:num_samples
            idx_vec = randperm(num_atoms, k);
            Z(idx_vec, i) = 1;
        end
        
    case 'bernoulli_gaussian'
        Z = rand(num_atoms, num_samples) <= (k/num_atoms);
end

A = randn(num_atoms, num_samples);
A(Z==0) = 0;

X = U*A;