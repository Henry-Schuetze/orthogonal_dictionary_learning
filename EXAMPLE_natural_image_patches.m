% clear;

% load gray value image to extract patches from
img = imread('test_images/barbara.png');
img = im2double(img);

% patch height and patch width
patch_size = 8;
num_dims = patch_size^2;

% the shift by which patches are extracted from the image
stride = 4;     % set stride = patch_size to get disjoint patches
                % set stride = 1 to get all overlapping patches

% CHOOSE A METHOD...
% ===============
%   either from CONSTRAINED K-SPARSE model:
%       {'CA', 'OSC', 'GF-OSC_cooling_learnrate', 'GF-OSC_line_search'}
%
%   or from UNCONSTRAINED LAGRANGIAN model:
%       {'DDTFC',  'lambda-OSC', 'lambda-GF-OSC_cooling_learnrate', 'lambda-GF-OSC_line_search'} 
method = 'OSC';
params = default_learning_params(method);
params.num_epochs = 200;

% Choose the user sparsity parameter depending on params.sparse_mode
params.sparsity_param = 8;      % constrained model
% params.sparsity_param = .2;     % unsonstrained lagrangian model

params.X = crop_patches_of_img(img, patch_size, stride);

% initialize with analytic orthogonal dictionary
% load('bases/DCT_08x08.mat', 'U_ref');
% U_init = U_ref;
% clear('U_ref');

% initialize with random orthogonal dictionary
U_init = solve_orth_procrustes(randn(num_dims));

params.U_init = U_init;

figure(1)
set(gcf, 'units', 'normalized', 'outerposition', [.5 .5 .5 .5])
clf;

subplot(1,3,1);
imshow(img);
title('image');

subplot(1,3,2);
result = orthogonal_dictionary_learning(params);
title('U_{final}');

subplot(1,3,3);
plot(result.cost_vec);
title('costs during epochs');
axis tight;
ylim([0 max(result.cost_vec)]);
xlabel('epoch');
ylabel('Cost Function Value')
grid on;