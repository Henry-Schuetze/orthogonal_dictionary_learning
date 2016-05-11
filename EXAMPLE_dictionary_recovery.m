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

% clear;

% Choose the reference basis from a combination of: 
% {'DCT', 'Haar_std', 'Haar_non-std'} x {'04x04', '08x08', '16x16'}
load('bases/DCT_08x08.mat', 'U_ref');

num_dims = size(U_ref,1);
img = create_dictionary_image(U_ref, [], true);

num_samples = 1000;
sparsity_param = 10;
synthetic_mode = 'column_k-sparse';

[X, Z] = create_synthetic_data_set(U_ref, num_samples, sparsity_param, synthetic_mode);

% initialize dictionary by a random orthogonal dictionary
num_atoms = num_dims;
U_init = solve_orth_procrustes(randn(num_dims, num_atoms));

%%
clear('params');

% CHOOSE A METHOD...
% ===============
%   either from CONSTRAINED K-SPARSE model:
%       {'CA', 'OSC', 'GF-OSC_cooling_learnrate', 'GF-OSC_line_search'}
%
%   or from UNCONSTRAINED LAGRANGIAN model:
%       {'DDTFC',  'lambda-OSC', 'lambda-GF-OSC_cooling_learnrate',
%       'lambda-GF-OSC_line_search'} 
method = 'GF-OSC_cooling_learnrate';
params = default_learning_params(method);

% ==> !IMPORTANT! <==: Choose sparsity parameter depending on params.sparse_mode
params.sparsity_param = sparsity_param;       % constrained model: expected number of non-zero coefficients per sample
% params.sparsity_param = .4;                  % unsonstrained model: hard threshold

params.X = X;
params.U_ref = U_ref;
params.U_init = U_init;

% uncomment the following lines if you wish to iteratively create image
% files of the dictionary
%
% params.write_dict_img_interval = 1;
% params.dict_img_path = '/home/some/path/'; % let the string end with '/'

figure(1)
set(gcf, 'units', 'normalized', 'outerposition', [.5 .5 .5 .5])
clf;

subplot(2,2,1);
imshow(img);
title('U_{ref} (ground truth)');

subplot(2,2,2);
result = orthogonal_dictionary_learning(params);
title('U_{final}');

subplot(2,2,3);
plot(result.cost_vec);
title('costs during epochs');
axis tight;
ylim([-0.05 max(result.cost_vec)]);
xlabel('epoch');
ylabel('Cost Function Value')
grid on;

subplot(2,2,4);
errorbar(mean(result.sim_mat,1), var(result.sim_mat,[],1));
title('similarity(U, U_{ref})');
axis tight;
ylim([-.05 1.05]);
xlabel('epoch');
ylabel('Mean Max Overlap');
grid on;