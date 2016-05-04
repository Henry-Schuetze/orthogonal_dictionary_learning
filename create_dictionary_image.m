% dict_img = create_dictionary_image(C, atom_seq_vec, discard_dc_flag, ...
%   grid_line_color, normalize_patch_flag)
%
% Creates an image of the (not necessarily orthogonal) dictionary C of
% which each atom (column) is visualized as a square patch arranged on a
% (square) grid. The number of rows of matrix C (dimensionality) has to be 
% a square number. By default, each patch is normalized for visualization,
% i.e., its mean value is subtracted, subsequently the patch is scaled to 
% unit supremum norm, and finally the resulting patch range [-1, 1] is 
% mapped to the gray value interval [0, 1] (by shifting and scaling).
%
% INPUT:
% ======
% C (required):
%   dictionary (num_dims x num_atoms)
%
% atom_seq_vec (optional, default: 1:num_atoms):
%   index sequence that defines the order of atom patches in the dictionary
%   image. atom_seq_vec should be a permutation of 1:num_atoms. The patches
%   will be arranged on the grid in a "column major" fashion.
%
% discard_dc_flag (optional, default: false):
%   if true, the atom of C that is closest to the dc component is
%   identified. The corresponding patch is only shifted by the average gray
%   value (0.5) and not further processed. Each column of C is assumed to
%   have unit length.
%
% grid_line_color (optional, default: 1):
%   defines the color of grid lines between the patches; possible values
%   are: 0 (black), 1 (red), 2 (green) or 3 (blue)
% 
% normalize_patch_flag (optional, default: true)
%   if true, patches are allowed to be normalized, i.e., modified for
%   visualization purpose (see function normalize_patch)
%
% OUTPUT:
% =======
% dict_img:
%   image of the dictionary. Each atom is represented by a suare patch. The
%   atom patches are ordered in a "column major" fashion according to 
%   atom_seq_vec.

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
function dict_img = create_dictionary_image(C, atom_seq_vec, ...
    discard_dc_flag, grid_line_color, normalize_patch_flag)

[num_dims, num_atoms] = size(C);

if nargin < 2 || isempty(atom_seq_vec)
    atom_seq_vec = 1:num_atoms;
end

if nargin < 3 || ~islogical(discard_dc_flag)
    discard_dc_flag = false;
end

if nargin < 4
    grid_line_color = 1;
end

if nargin < 5
    normalize_patch_flag = true;
end

if discard_dc_flag
    % estimate which atom is most similar to dc component
    dc_vec = (1/sqrt(num_dims))*ones(num_dims, 1);
    [~, dc_idx] = max(abs(C(:,atom_seq_vec)' * dc_vec));
    [dc_x, dc_y] = ind2sub(ceil(sqrt(num_atoms))*[1, 1], dc_idx);
end

if grid_line_color ~= 0
    dict_img = ones([ceil(sqrt(num_atoms))*(1+sqrt(num_dims))-1, ceil(sqrt(num_atoms))*(1+sqrt(num_dims))-1, 3]);
    
    switch grid_line_color
        case 1  % red
            color_channel = [2 3];
        case 2  % green
            color_channel = [1 3];
        case 3  % blue
            color_channel = [1 2];
        otherwise
            error('unexpected value of grid_line_color');
    end
    dict_img(:,:,color_channel) = 0;
else
    dict_img = nan(ceil(sqrt(num_atoms))*(1+sqrt(num_dims))-1);
end

for atom_idx = 1:ceil(sqrt(num_atoms))^2

    [x, y] = ind2sub(ceil(sqrt(num_atoms))*[1, 1], atom_idx);
    
    if atom_idx <= length(atom_seq_vec)
        atom = C(:, atom_seq_vec(atom_idx));

        patch = reshape(atom, sqrt(num_dims), sqrt(num_dims));
        if discard_dc_flag && (x==dc_x) && (y==dc_y)
            % shift dc component by mean gray value
            patch = patch + 0.5;
        else
            if normalize_patch_flag
                patch = normalize_patch(patch);
            end
        end
    else
        patch = zeros(sqrt(num_dims));
    end

    if grid_line_color ~= 0
        % colored grid lines
        patch = repmat(patch, [1 1 3]);
        dict_img(((x-1)*(sqrt(num_dims)+1)+1) : x*(sqrt(num_dims)+1)-1, ...
            ((y-1)*(sqrt(num_dims)+1)+1) : y*(sqrt(num_dims)+1)-1, :) = patch;
    else
        % black grid lines
        dict_img(((x-1)*(sqrt(num_dims)+1)+1) : x*(sqrt(num_dims)+1)-1, ...
            ((y-1)*(sqrt(num_dims)+1)+1) : y*(sqrt(num_dims)+1)-1) = patch;
    end
end

function patch = normalize_patch(patch)

patch = patch - mean(patch(:));
patch = patch ./ max(abs(patch(:)));
patch = (patch + 1) / 2;