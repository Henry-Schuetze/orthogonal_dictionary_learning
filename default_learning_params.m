% params = default_learning_params(method)
%
% Sets default parameters to learn an orthogonal sparse coding dictionary 
% by a particular method. The function returns a struct with fields: for
% instance a callback function .update() or parameters depending on the
% particular method. The parameter struct is intended to be extended (e.g.
% by data matrix X), and subsequently passed to the main learning function
% orthogonal_dictionary_learning.m. Remember that caling 
% default_learning_params(method) sets a 'fresh' seed for the random number
% generator (params.rand_seed). If you intend a deterministic random 
% initialization of the dictionary and a deterministic sequence of training
% data samples to be processed (the latter only for online methods) by the
% main learning function, overwrite the field .rand_seed.
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
function params = default_learning_params(method)

params = struct( ...
	'num_epochs', 1e2, ...          % number of learning epochs
    'sim_interval', 1, ...          % interval to compute similarity to reference dictionary
    'cost_interval', 1, ...         % interval to evaluate cost function
    'show_dict_img_interval', 1, ...% interval to visualize dictionary
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
        % use update_OSC_fast(params) instead of update_OSC(params)
        % params.update = @(params) update_OSC(params);
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
        params.tau = 5e-1;
        params.c = 5e-1;
        params.alpha0 = 5e-2;
        params.update = @(params) update_GFOSC_line_search(params);
               
    otherwise
        error('unknown method: %s', method);
end

params.rand_seed = cputime;