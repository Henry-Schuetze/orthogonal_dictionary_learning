% params = default_learning_params(method)
%
% Sets default parameters to learn an orthogonal sparse coding dictionary 
% by a particular method. The function returns a struct which is intended 
% to be extended (e.g. by data), and subsequently passed to the main
% learning function orthogonal_dictionary_learning.m.
% 
% INPUT:
% ======
% method (required):
%   a string selecting the learning method
% 
% OUTPUT:
% =======
% params
%   a struct containing the parameter values

% Henry Schuetze
% Institute for Neuro- and Bioinformatics
% University of Luebeck, Germany
% Henry.Schuetze@uni-luebeck.de
function params = default_learning_params(method)

params = struct( ...
	'num_epochs', 1e2, ...          % number of learning epochs
    'sim_interval', 1, ...          % interval to compute similarity to reference dictionary
    'cost_interval', 1, ...         % interval to evaluate cost function
    'plot_dict_interval', 1, ...    % interval to visualize dictionary
    'sim_stop_thresh', 0.9999, ...  % if similarity exceed this threshold, learning is stopped
    'method', method, ...           % learning method
    'verbose_flag', true ...        % to print messages during learning
);

switch method
	case 'CA'
        params.learn_type = 'batch';
        params.sparse_mode = 'column_k-sparse';
        params.update = @(params) update_orth_procrustes(params);
        
    case 'DDTFC'
        params.learn_type = 'batch';
        params.sparse_mode = 'hard_thresh';
        params.update = @(params) update_orth_procrustes(params);

    case 'OSC'
        params.learn_type = 'online';
        params.sparse_mode = 'column_k-sparse';
        params.eps_i = 1e-1;
        params.eps_f = 1e-3;
        params.update = @(params) update_OSC_fast(params);
        
    case 'lambda-OSC'
        params.learn_type = 'online';
        params.sparse_mode = 'hard_thresh';
        params.eps_i = 1e-1;
        params.eps_f = 1e-3;
        params.update = @(params) update_lambdaOSC_fast(params);

    case 'GF-OSC_cooling_learnrate'
        params.learn_type = 'online';
        params.sparse_mode = 'column_k-sparse';
        params.eps_i = 5e-2;
        params.eps_f = 5e-3;
        params.update = @(params) update_GFOSC_cooling_learnrate(params);
        
    case 'GF-OSC_line_search'
        params.learn_type = 'online';
        params.sparse_mode = 'column_k-sparse';
        params.tau = .5;
        params.c = .5;
        params.alpha0 = 2;
        params.update = @(params) update_GFOSC_line_search(params);
        
    case 'lambda-GF-OSC_cooling_learnrate'
        params.learn_type = 'online';
        params.sparse_mode = 'hard_thresh';
        params.eps_i = 5e-2;
        params.eps_f = 5e-3;
        params.update = @(params) update_GFOSC_cooling_learnrate(params);
        
    case 'lambda-GF-OSC_line_search'
        params.learn_type = 'online';
        params.sparse_mode = 'hard_thresh';
        params.tau = .5;
        params.c = .5;
        params.alpha0 = 5e-2;
        params.update = @(params) update_GFOSC_line_search(params);
               
    otherwise
        error('unknown method: %s', method);
end

if strcmp(params.learn_type, 'online')
    params.rand_seed = cputime;
end