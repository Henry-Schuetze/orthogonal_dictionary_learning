% [sim_vec, match_atoms_vec] = dictionary_similarity(C, C_ref)
%
% Computes the overlaps between "matching" columns of matrix C and matrix 
% C_ref. The matching is assigned based on 1st, 2nd , 3rd, ... largest
% overlap. The assignment is a bijective mapping. The columns of C and
% C_ref should all have unit length.
%
% INPUT:
% ======
% C (required):
%   dictionary (num_dims x num_atoms)
%
% C_ref (required):
%   reference dictionary (num_dims x num_atoms)
%
% OUTPUT:
% =======
% sim_vec:
%   vector (1 x num_atoms) containing the overlaps of matched atom pairs
%
% match_atoms_vec:
%   defines column permutation of matrix C such that a column
%   C(:, match_atoms_vec(col_idx)) matches with column C_ref(:, col_idx)
%

% Henry Schuetze 
% Institute for Neuro- and Bioinformatics
% University of Luebeck, Germany
% Henry.Schuetze@uni-luebeck.de
function [sim_vec, match_atoms_vec] = dictionary_similarity(C, C_ref)

num_dims = size(C, 1);
G = abs(C_ref'*C);

match_atoms_vec = zeros(1, num_dims);
sim_vec = zeros(1, num_dims);

for h = 1:num_dims
    [~, idx] = max(G(:));
    [i,j] = ind2sub(size(G), idx);
    
    match_atoms_vec(i) = j;
    sim_vec(h) = G(i, j);

    G(i,:) = 0;
    G(:,j) = 0;
end