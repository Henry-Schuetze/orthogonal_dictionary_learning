clear;

% Choose the reference basis from a combination of: 
% {'DCT', 'Haar_std', 'Haar_non-std'} x {'04x04', '08x08', '16x16'}
load('bases/DCT_08x08.mat', 'U_ref');

num_dims = size(U_ref,1);
img = create_dictionary_image(U_ref, [], true);

num_samples = 1000;
sparsity_param = 10;
synthetic_mode = 'column_k-sparse';

[X, Z] = create_synthetic_data_set(U_ref, num_samples, sparsity_param, synthetic_mode);

% create an initial random basis
U_init = solve_orth_procrustes(randn(num_dims));

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

% <IMPORTANT>: Choose sparsity parameter depending on params.sparse_mode
params.sparsity_param = sparsity_param;       % constrained model: expected number of non-zero coefficients per sample
% params.sparsity_param = 0.4;                    % unsonstrained model: hard threshold

params.X = X;
params.U_ref = U_ref;
params.U_init = U_init;

% uncomment the following lines if you wish to iteratively create image
% files of the dictionary 
%
% params.img_seq_interval = 1;
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
ylim([0 max(result.cost_vec)]);
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