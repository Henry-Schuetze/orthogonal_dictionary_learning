% dict_img = create_dictionary_image(C, seq_vec, discard_dc_flag, ...
%   grid_line_color, scale_flag)
%
% Creates a image of dictionary C in which each atom is represented as a 
% square patch arranged on a (square) grid. The number of rows of matrix C 
% has to be a square number. From each patch its mean is subtracted, 
% subsequently the patch is scaled to unit supremum norm and finally 
% shifted and rescaled to match the interval [0 1].
%
% INPUT:
% ======
% C (required):
%   dictionary (num_dims x num_atoms)
%
% seq_vec (optional, default: 1:num_atoms):
%   index sequence that defines the order of atoms in the dictionary image.
%   seq_vec should be a permutation of 1:num_atoms. The patches will be
%   arranged on the grid in a "column major" fashion.
%
% discard_dc_flag (optional, default: false):
%   if true the dc component is estimated and its scaling to unit supremum 
%   norm is prevented
%
% grid_line_color (optional, default: 1):
%   defines the color of grid lines between the patches; possible values
%   are: 0 (black), 1 (red), 2 (green) or 3 (blue)
% 
% scale_flag (optional, default: true)
%   if true patches are allowed to be normalized (see normalize_patch)
%
% OUTPUT:
% =======
% dict_img:
%   image of the dictionary. Each atom is represented by a suare patch. The
%   atom patches are ordered in column major according to seq_vec.

% Henry Schuetze 
% Institute for Neuro- and Bioinformatics
% University of Luebeck, Germany
% Henry.Schuetze@uni-luebeck.de
function dict_img = create_dictionary_image(C, seq_vec, ...
    discard_dc_flag, grid_line_color, scale_flag)

[num_dims, num_atoms] = size(C);

if nargin < 2 || isempty(seq_vec)
    seq_vec = 1:num_atoms;
end

if nargin < 3 || ~islogical(discard_dc_flag)
    discard_dc_flag = false;
end

if nargin < 4
    grid_line_color = 1;
end

if nargin < 5
    scale_flag = true;
end

if discard_dc_flag
    % estimate which atom is most similar to dc component
    dc_vec = (1/sqrt(num_dims))*ones(num_dims, 1);
    [~, dc_idx] = max(abs(C(:,seq_vec)' * dc_vec));
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
    
    if atom_idx <= length(seq_vec)
        atom = C(:, seq_vec(atom_idx));

        if discard_dc_flag && (x==dc_x) && (y==dc_y)
            % shift dc component by mean gray value
            patch = reshape(atom, sqrt(num_dims), sqrt(num_dims)) + 0.5;
        else
            patch = reshape(atom, sqrt(num_dims), sqrt(num_dims));
            if scale_flag
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