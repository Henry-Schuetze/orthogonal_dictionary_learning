function U = solve_orth_procrustes(M)

[num_dims, num_atoms] = size(M);
assert(num_atoms <= num_dims);

[Q, S, R] = svd(M);

% % set all singular values to one
S = eye(size(S));

U = Q*S*R';