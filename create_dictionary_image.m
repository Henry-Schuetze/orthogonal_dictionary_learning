% dict_img = create_dictionary_image(C, seq_vec, discard_dc_flag, ...
%   grid_line_color, scale_flag)
%
% Creates a image of dictionary C in which each atom is represented as a 
% square patch on a (square) grid. Thus the number of rows of matrix C 
% needs to be a square number. From each patch its mean is subtracted,
% subsequently the patch is scaled to unit supremum norm and finally 
% shifted and rescaled into the interval [0 1].
%
% INPUT:
% ======
% C (required):
%   dictionary (num_dims x num_atoms)
%
% seq_vec (optional, default: 1:num_atoms):
%   sequence that defines the order of atoms in the dictionary image. 
%   seq_vec should be a permutation of 1:num_atoms.
%
% discard_dc_flag (optional, default: false):
%   boolean to estimate the dc component and prevent its scaling to unit
%   supremum norm
%
% grid_line_color (optional, default: 1):
%   defines the color of grid lines between the patches in black
%   (0), red (1), green (2) or blue (3)
% 
% scale_flag (optional, default: true)
%   defines if patches are allowed to be scaled (i.e. normalized)
%
% OUTPUT:
% =======
% dict_img:
%   image of the dictionary. Each atom is represented by a suare patch. The
%   atom patches are ordered in column major according to seq_vec.
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

for i = 1:ceil(sqrt(num_atoms))^2

    [x, y] = ind2sub(ceil(sqrt(num_atoms))*[1, 1], i);
    
    if i <= length(seq_vec)
        u = C(:, seq_vec(i));

        if discard_dc_flag && (x==dc_x) && (y==dc_y)
            patch = reshape(u, sqrt(num_dims), sqrt(num_dims)) + 0.5;
        else
            patch = reshape(u, sqrt(num_dims), sqrt(num_dims));
            if scale_flag
                patch = normalize_patch(patch);
            end
        end
    else
        patch = zeros(sqrt(num_dims));
    end

    if grid_line_color ~= 0
        patch = repmat(patch, [1 1 3]);
        dict_img(((x-1)*(sqrt(num_dims)+1)+1) : x*(sqrt(num_dims)+1)-1, ...
            ((y-1)*(sqrt(num_dims)+1)+1) : y*(sqrt(num_dims)+1)-1, :) = patch;
    else
        dict_img(((x-1)*(sqrt(num_dims)+1)+1) : x*(sqrt(num_dims)+1)-1, ...
            ((y-1)*(sqrt(num_dims)+1)+1) : y*(sqrt(num_dims)+1)-1) = patch;
    end
end

function patch = normalize_patch(patch)

patch = patch - mean(patch(:));
patch = patch ./ max(abs(patch(:)));
patch = (patch + 1) / 2;