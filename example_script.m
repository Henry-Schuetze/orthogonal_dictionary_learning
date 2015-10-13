clear;

% Choose the reference basis from a combination of: 
% {'DCT', 'Haar_std', 'Haar_non-std'} x {'04x04', '08x08', '16x16'}
load('bases/DCT_08x08.mat', 'U_ref');

num_dims = size(U_ref,1);
img = create_dictionary_image(U_ref, [], true);

num_samples = 1000;
k = 12;
synth_mode = 'column_k-sparse';

[X, Z] = create_synthetic_data_set(U_ref, num_samples, k, synth_mode);

U_init = solve_orth_procrustes(randn(num_dims));

%%
clear('params');

% Choose a method...
%
%   either from CONSTRAINED K-SPARSE MODEL:
%       {'CA', 'OSC', 'GF-OSC_cooling_learnrate', 'GF-OSC_line_search'}
%
%   or from UNCONSTRAINED LAGRANGIAN MODEL:
%       { 'DDTFC',  'lambda-OSC', 'lambda-GF-OSC_cooling_learnrate', 'lambda-GF-OSC_line_search'} 
method = 'lambda-GF-OSC_cooling_learnrate';
params = default_learning_params(method);

% Choose the user sparsity parameter depending on params.sparse_mode
% params.sparsity_param = k;      % constrained model
params.sparsity_param = .8;     % unsonstrained lagrangian model

params.X = X;
params.U_ref = U_ref;
params.U_init = U_init;

figure(1)
set(gcf,'units','normalized','outerposition',[.5 .5 .5 .5])
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