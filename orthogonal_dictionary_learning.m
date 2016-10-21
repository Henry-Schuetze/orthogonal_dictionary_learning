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
% params.rand_seed (optional, default: 0):
%   a seed for the random number generator such that the randomly generated
%   initial dictionary and the sequence of training data samples for online
%   methods are deterministic
%
% params.U_init (optional):
%   initial orthogonal dictionary. Provide params.U_init with less columns
%   than rows in order to learn an undercomplete dictionary
%
% params.cost_interval (optional):
%   if set, a cost function is evaluated with respect to X and U each time
%   the epoch counter has a value that is a multiple of cost_interval
%
% params.U_ref (optional):
%   an orthogonal reference dictionary, required to compute the similarity
%   between iterating dictionary and a ground-truth dictionary
%
% params.sim_interval (optional, default: 1):
%   iff U_ref is set, the similarity between U and U_ref is computed each 
%   time the epoch counter has a value that is a multiple of sim_interval
%
% params.sim_stop_thresh (optional, default: .9999):
%   iff U_ref is set, the iteration is stopped as soon as the epoch counter
%   has a value that is a multiple of sim_interval and the smallest overlap
%   among all matching pairs of atoms from U and U_ref is above
%   sim_stop_thresh.
%
% params.write_dict_img_interval (optional):
%   if set, an image file of the current dictionary is written (assuming
%   the atoms represent square patches) each time the epoch counter has a
%   value that is a multiple of write_dict_img_interval. A path can be
%   supplied via the string params.dict_img_path, otherwise the images are
%   stored into the current path
%
% params.show_dict_img_interval (optional, default: 1): 
%   if set, an image of the current dictionary is shown (assuming the
%   atoms represent square patches) each time the epoch counter has a value
%   that is a multiple of show_dict_img_interval.
%
% params.verbose_flag (optional, default: true)
%   if set, output messages are printed during iterations
%
% OUTPUT: result
% =======
% result.U:
%   final orthogonal dictionary
%
% result.cost_vec:
%   cost_function values for every params.cost_interval epoch
%
% result.sim_mat:
%   matrix that columnwise contains the overlaps of matching atoms from U
%   and U_ref

% Copyright © 2016 Henry Schuetze
% Institute for Neuro- and Bioinformatics
% University of Luebeck, Germany
% Henry.Schuetze@uni-luebeck.de
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
function result = orthogonal_dictionary_learning(params)

if ~isfield(params, 'rand_seed')
    params.rand_seed = 0;
end
% set seed of the random number generator
rng(params.rand_seed, 'twister');

[num_dims, num_samples] = size(params.X);

if ~isfield(params, 'U_init')
    % randomly create initial orthogonal dictionary
    params.U_init = solve_orth_procrustes(randn(num_dims));
end
params.U = params.U_init;

num_atoms = size(params.U, 2);

% prepare iterative computation of cost function
cost_flag = isfield(params, 'cost_interval');
if cost_flag
    params.cost_vec = zeros(1, floor(params.num_epochs/params.sim_interval)+1);
end

% prepare iterative computation of similarity between U and U_ref
sim_flag = isfield(params, 'U_ref');
if sim_flag
    if ~isfield(params, 'sim_interval')
        params.sim_interval = 1;
    end
    
    if ~isfield(params, 'sim_stop_thresh')
        params.sim_stop_thresh = .9999;
    end

    params.sim_mat = zeros(num_atoms, floor(params.num_epochs/params.sim_interval)+1);
end

% prepare iterative creation of dictionary image files
write_dict_img_flag = isfield(params, 'write_dict_img_interval');
if write_dict_img_flag
    if ~isfield(params, 'dict_img_path')
        params.dict_img_path = '';
    end
    
    dict_img = create_dictionary_image(params.U);
    dict_img_filename = sprintf('dictionary_image_epoch_%0.8d.png', 0);
    imwrite(dict_img, [params.dict_img_path, dict_img_filename], 'png');
end

show_dict_img_flag = isfield(params, 'show_dict_img_interval');

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
end

for epoch = 0:params.num_epochs
    
    % do not update in epoch 0, only visualize and store subject to initial
    % dictionary
    if epoch > 0
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
    end
    
    if (show_dict_img_flag && ~mod(epoch, params.show_dict_img_interval)) ...
            || (write_dict_img_flag && ~mod(epoch, params.write_dict_img_interval))
        % if U_ref exists
        if sim_flag
            % call function dictionary_similarity to match the atoms of U
            % with the atoms of U_ref for the dictionary image
            [~, match_atoms_vec] = dictionary_similarity(params.U, params.U_ref);
        else
            match_atoms_vec = [];
        end
        
        % create dictionary image
        dict_img = create_dictionary_image(params.U, match_atoms_vec, true);
        
        if show_dict_img_flag && ~mod(epoch, params.show_dict_img_interval)
            % show dictionary image
            imshow(dict_img);
            title(sprintf('U after epoch %08d', epoch));
            drawnow;
        end
        
        if write_dict_img_flag && ~mod(epoch, params.write_dict_img_interval)
            % write dictionary image to file system
            dict_img_filename = sprintf('dictionary_image_epoch_%0.8d.png', epoch);
            imwrite(dict_img, [params.dict_img_path, dict_img_filename], 'png');
        end
    end
    
    % evaluate cost_function
    if cost_flag && ~mod(epoch, params.cost_interval)
        params.cost_vec((epoch/params.cost_interval)+1) = ...
            cost_function(params.X, params.U, params.sparse_mode, params.sparsity_param);
    end
    
    % compute similarity between params.U and params.U_ref
    if sim_flag && ~mod(epoch, params.sim_interval)
        params.sim_mat(:,(epoch/params.sim_interval)+1) = ...
            dictionary_similarity(params.U, params.U_ref);
        
        % if similarity stop condition is satisfied, resize params.sim_mat
        % and params.cost_vec
        if min(params.sim_mat(:, floor(epoch/params.sim_interval)+1)) >= params.sim_stop_thresh
            params.sim_mat = params.sim_mat(:, 1:(floor(epoch/params.sim_interval)+1));
            if cost_flag
                params.cost_vec = params.cost_vec(1:(floor(epoch/params.cost_interval)+1));
            end
            
            if verbose_flag
                fprintf('similarity stop condition reached. finished after %d epochs.\n', epoch);
            end
            break;
        end
    end
    
    if verbose_flag && epoch > 0
        fprintf('epoch %08d completed.\n', epoch);
    end
end

if ~isfield(params, 'A')
    params.A = sparse_coefficients(params.X, params.U, params.sparse_mode, params.sparsity_param);
end

result = params;

if isfield(result, 'show_dict_img_interval')
    result = rmfield(result, 'show_dict_img_interval');
end
if isfield(result, 'verbose_flag')
    result = rmfield(result, 'verbose_flag');
end
if isfield(result, 'x')
    result = rmfield(result, 'x');
end
if isfield(result, 't')
    result = rmfield(result, 't');
end
if isfield(result, 't_max')
    result = rmfield(result, 't_max');
end