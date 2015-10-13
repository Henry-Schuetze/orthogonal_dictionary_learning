function params = update_orth_procrustes(params)

params.A = sparse_coefficients(params.X, params.U, ...
    params.sparse_mode, params.sparsity_param);
params.U = solve_orth_procrustes(params.X * params.A');