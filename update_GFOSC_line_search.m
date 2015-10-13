function params = update_GFOSC_line_search(params)

alpha = params.alpha0;

a_hat = sparse_coefficients(params.x, params.U, params.sparse_mode, params.sparsity_param);

x_hat = params.U*a_hat;
grad_E_A = x_hat*params.x' - params.x*x_hat';

t = -params.c*(-grad_E_A(:)' * grad_E_A(:));
cost_at_U = cost_function(params.x, params.U, params.sparse_mode, params.sparsity_param);

delta_U = expm(alpha*(-grad_E_A));
cost_at_delta_UU = cost_function(params.x, delta_U*params.U, params.sparse_mode, params.sparsity_param);

while (cost_at_U - cost_at_delta_UU) < alpha*t
    alpha = params.tau*alpha;
    delta_U = expm(-alpha*grad_E_A);
    cost_at_delta_UU = cost_function(params.x, delta_U*params.U, ...
        params.sparse_mode, params.sparsity_param);
end

params.U = delta_U*params.U;
