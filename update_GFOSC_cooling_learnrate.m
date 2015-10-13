function params = update_GFOSC_cooling_learnrate(params)

eps_t = params.eps_i * (params.eps_f/params.eps_i)^(params.t/params.t_max);

a_hat = sparse_coefficients(params.x, params.U, params.sparse_mode, ...
    params.sparsity_param);
x_hat = params.U*a_hat;

grad_E_A = x_hat*params.x' - params.x*x_hat';
delta_U = expm(-eps_t*grad_E_A);

params.U = delta_U*params.U;