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

clear;

% load gray value image to extract patches from
% img = imread('test_images/boat.png');
% img = imread('test_images/pirate.tif');
img = imread('test_images/barbara.png');
img = im2double(img);

% CHOOSE A METHOD...
% ===============
%   either from CONSTRAINED K-SPARSE model:
%       {'CA', 'OSC', 'GF-OSC_cooling_learnrate', 'GF-OSC_line_search'}
%
%   or from UNCONSTRAINED LAGRANGIAN model:
%       {'DDTFC',  'lambda-OSC', 'lambda-GF-OSC_cooling_learnrate', 
%       'lambda-GF-OSC_line_search'} 
method = 'CA';
params = default_learning_params(method);

% ==> !IMPORTANT! <==: choose sparsity parameter depending on params.sparse_mode
params.sparsity_param = 10;      % constrained model: expected number of non-zero coefficients per sample
% params.sparsity_param = .2;     % unsonstrained model: hard threshold

params.num_epochs = 210;
% params.rand_seed = 0;

% set height (= width) of patches to extract from image
patch_size = 8;
num_dims = patch_size^2;

% set shift/stride by which patches are extracted from the image
stride = 8;     % set stride = patch_size to get disjoint patches
                % set stride = 1 to get all overlapping patches (increases
                % number of samples to approximately the number of pixels)
params.X = crop_patches_of_img(img, patch_size, stride);

% initialize dictionary by an analytic orthogonal dictionary
% load('bases/DCT_08x08.mat', 'U_ref');
% U_init = U_ref;
% clear('U_ref');

% initialize dictionary by a random orthogonal dictionary
num_atoms = num_dims;
U_init = solve_orth_procrustes(randn(num_dims, num_atoms));
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

% semilog plot of costs including initial conditions (epoch 0)
subplot(1,3,3);
semilogx(result.cost_vec);
h = gca;
minor_xtick_vec = 0;
for i = 1:ceil(log10(result.cost_vec))
    minor_xtick_vec = [minor_xtick_vec, (1:10) .* 10^(i-1)];
end
minor_xtick_vec =  unique(minor_xtick_vec+1);
h.XAxis.MinorTickValues = unique(minor_xtick_vec);
major_xtick_vec = [0, 10.^(0:ceil(log10(result.cost_vec)))]+1;
set(gca, 'XTick', major_xtick_vec, 'XtickLabel', major_xtick_vec-1);
title('costs during epochs');
axis tight;
ylim([0 max(result.cost_vec)]);
xlabel('epoch');
ylabel('Cost Function Value')
grid on;