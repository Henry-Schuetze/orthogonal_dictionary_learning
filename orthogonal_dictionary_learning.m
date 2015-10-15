% function result = orthogonal_dictionary_learning(params)
%
% INPUT: params
% ======
% params.X (required):
%   data matrix (num_dims x num_samples), where each column represents a
%   sample
%
% params.sparse_mode (required):
%   a string with value either 'column_k-sparse' or 'hard_thresh' selecting
%   the sparse model
%
% params.sparsity_param (required):
%   sparsity parameter corresponding to the sparse model, i.e., either the 
%   number of non-zero coefficients (if sparse_mode == 'column_k-sparse') 
%   or the hard threshold (if sparse_mode == 'hard_thresh')
%
% params.num_epochs (required):
%   number of full batch iterations
%
% params.U_init (optional):
%   initial complete orthogonal dictionary
%
% params.cost_interval (optional):
%   if set, a cost function is evaluated with respect to X and U each time
%   the epoch counter has a value which is a multiple of cost_interval
%
% params.U_ref (optional):
%   a complete orthogonal reference dictionary, required to compute
%   similarity between iterating dictionary and a ground-truth dictionary
%
% params.sim_interval (optional, default: 1):
%   iff U_ref is set, the similarity between U and U_ref is computed each 
%   time the epoch counter has a value which is a multiple of cost_interval
%
% params.sim_stop_thresh (optional, default: .9999):
%   iff U_ref is set, the iteration is stopped as soon as the smallest 
%   overlap among all matching atom pairs from U and U_ref is above 
%   sim_stop_thresh
%
% params.img_seq_interval (optional):
%   if set, an image file of the current dictionary is written (assuming
%   the atoms represent square patches) each time the epoch counter has a
%   value which is a multiple of img_seq_interval
%
% params.plot_dict_interval (optional, default: false): 
%   if set, an image of the current dictionary is shown (assuming the
%   atoms represent square patches) each time the epoch counter has a value
%   which is a multiple of plot_dict_interval
%
% params.verbose_flag (optional, default: true)
%   if set, output messages are printed during iterations
%
% OUTPUT: result
% =======
% result.U:
%   final complete orthogonal dictionary
%
% result.cost_vec:
%   cost_function values for every params.cost_interval epoch
%
% result.sim_mat:
%   matrix that columnwise contains the overlaps of matching atoms from U
%   and U_ref

% Henry Schuetze 
% Institute for Neuro- and Bioinformatics
% University of Luebeck, Germany
% Henry.Schuetze@uni-luebeck.de
function result = orthogonal_dictionary_learning(params)

result = struct('U', []);

[num_dims, num_samples] = size(params.X);

if isfield(params, 'U_init')
    % load initial dictionary if supplied
    params.U = params.U_init;
else
    % randomly create initial orthogonal dictionary
    params.U = solve_orth_procrustes(randn(num_dims));
end

% prepare iterative computation of cost function
cost_flag = isfield(params, 'cost_interval');
if cost_flag
    result.cost_vec = zeros(1, floor(params.num_epochs/params.sim_interval));
end

% prepare iterative computation of similarity between U and U_ref
sim_flag = isfield(params, 'U_ref');
if sim_flag
	if ~isfield(params, 'sim_interval')
        params.sim_interval = 1;
    end

    result.sim_mat = zeros(num_dims, floor(params.num_epochs/params.sim_interval)+1);
	result.sim_mat(:,1) = dictionary_similarity(params.U_ref, params.U);
        
    if ~isfield(params, 'sim_stop_thresh')
        params.sim_stop_thresh = .9999;
    end
end

% prepare iterative creation of dictionary images
img_seq_flag = isfield(params, 'img_seq_interval');
if img_seq_flag
    dict_img = create_dictionary_image(params.U);
    filename = sprintf('dictionary_image_epoch_%0.8d.png', 0);
    imwrite(dict_img, filename, 'png');
end

plot_dict_flag = isfield(params, 'plot_dict_interval');

% set verbose_flag
if ~isfield(params, 'verbose_flag')
    verbose_flag = true;
else
    verbose_flag = params.verbose_flag;
end

% extend params for online learning methods
if strcmp(params.learn_type, 'online')
    % provide a learn step counter and the total number of learning steps 
    % in order to allow computation of variable parameters, e.g., a cooling
    % learning rate
    params.t = 0;
    params.t_max = params.num_epochs * num_samples;
    
    % set seed of the random number generator, e.g., to obtain a
    % deterministic random sequence of data samples
    if isfield(params, 'rand_seed')
        rng(params.rand_seed, 'twister');
    else
        rng(0, 'twister');
    end
end

for epoch = 1:params.num_epochs
    
    switch params.learn_type
        case 'batch'
            params = params.update(params);
            
        case 'online'
            for sample_idx = randperm(num_samples);
                params.x = params.X(:, sample_idx);
                params = params.update(params);
                params.t = params.t + 1;
            end
            
        otherwise
            error('unknown params.learn_type');
    end
    
    if (plot_dict_flag && ~mod(epoch, params.plot_dict_interval)) ...
            || (img_seq_flag && ~mod(epoch, params.img_seq_interval))
        
        if sim_flag
            % get sequence to match dictionry atoms with reference
            % dictionary
            [~, comp_seq_vec] = dictionary_similarity(params.U_ref, params.U);
        else
            comp_seq_vec = [];
        end
        
        % create dictionary image
        dict_img = create_dictionary_image(params.U, comp_seq_vec, true);
        
        if plot_dict_flag && ~mod(epoch, params.plot_dict_interval)
            imshow(dict_img);
            title(sprintf('U after epoch %08d', epoch));
            pause(.2);
        end
        
        % write dictionary image to file system
        if img_seq_flag && ~mod(epoch, params.img_seq_interval)
            filename = sprintf('dictionary_image_epoch_%0.8d.png', epoch);
            imwrite(dict_img, filename, 'png');
        end
    end
    
    % evaluate cost_function
    if cost_flag && ~mod(epoch, params.cost_interval)
        result.cost_vec(epoch/params.cost_interval) = ...
            cost_function(params.X, params.U, params.sparse_mode, params.sparsity_param);
    end
    
    % compute similarity between params.U_ref and params.U
    if sim_flag && ~mod(epoch, params.sim_interval)
        result.sim_mat(:, floor(epoch/params.sim_interval)+1) = ...
            dictionary_similarity(params.U_ref, params.U);
        
        % if stop condition is satisfied, resize result.sim_mat and
        % result.cost_vec
        if min(result.sim_mat(:, floor(epoch/params.sim_interval)+1)) >= params.sim_stop_thresh
            result.sim_mat = result.sim_mat(:, 1:(floor(epoch/params.sim_interval)+1));
            if cost_flag
                result.cost_vec = result.cost_vec(1:(floor(epoch/params.cost_interval)));
            end
            
            if verbose_flag
                fprintf('similarity stop condition reached. finished after %d epochs.\n', epoch);
            end
            break;
        end
    end
    
    if verbose_flag
        fprintf('epoch %08d completed.\n', epoch);
    end
end

result.U = params.U;
if ~isfield(params, 'A')
    params.A = sparse_coefficients(params.X, params.U, params.sparse_mode, params.sparsity_param);
end
result.A = params.A;