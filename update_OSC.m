function params = update_OSC(params)

[~, seq_vec] = sparse_coefficients(params.x, params.U, 'column_k-sparse', params);
eps_t = params.eps_i * (params.eps_f/params.eps_i)^(params.t/params.t_max);

x_r = params.x;

% tStart = tic;
for i = 1 : size(params.U,2)
    for l = 1 : (i-1)
        params.U(:, seq_vec(i)) = params.U(:, seq_vec(i)) - ...
            (params.U(:, seq_vec(i))'*params.U(:, seq_vec(l))) * params.U(:,seq_vec(l));
    end
    if i <= params.k
        % u_k = u_k / sqrt(sum(u_k.^2));
        params.U(:, seq_vec(i)) = params.U(:, seq_vec(i)) + eps_t * (params.U(:, seq_vec(i))'* x_r) * x_r;
    end
    params.U(:, seq_vec(i)) = params.U(:, seq_vec(i)) / sqrt(sum(params.U(:, seq_vec(i)).^2));
    x_r = x_r - (params.U(:, seq_vec(i))'*x_r)*params.U(:, seq_vec(i));
end

% tElapsed = toc(tStart);
% fprintf('dictionary update took: %.4fs\n', tElapsed);


% 
% function params.U = dictionary_update_OSC_II(params.U, x_r, seq_vec, K, eps_t)
% 
% tStart = tic;
% U_old = zeros(size(params.U,1), size(params.U,2)+1);
% U_old(:,1:(end-1)) = params.U;
% U_old(:,end) = x_r;
% 
% for i = 1 : size(params.U,2)
%     params.U(:, seq_vec(i)) = U_old(:, seq_vec(i)) + eps_t * (U_old(:, seq_vec(i))'* U_old(:,end));
%     params.U(:, seq_vec(i)) = params.U(:, seq_vec(i)) ./ sqrt(sum(params.U(:, seq_vec(i)).^2));
% 
%     U_old = U_old - params.U(:, seq_vec(i)) * (params.U(:, seq_vec(i))' * U_old);
%     
% %     params.U(:, seq_vec((i+1):end)) = params.U(:, seq_vec((i+1):end)) - ...
% %         params.U(:, seq_vec(i)) * (params.U(:, seq_vec(i))' *  params.U(:, seq_vec((i+1):end)));
%     
% %    x_r = x_r - ...
% %        U_new(:, seq_vec(i)) * (U_new(:, seq_vec(i))' *  x_r);
% end
% tElapsed = toc(tStart);
% fprintf('dictionary update took: %.4fs\n', tElapsed);