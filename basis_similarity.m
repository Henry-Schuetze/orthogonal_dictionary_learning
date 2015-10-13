% function [sim_vec, assign_vec] = basis_similarity(U_ref, U)
%
% Computes the overlaps between "matching" columns of matrix U and matrix 
% U_ref. The matching is assigned based on 1st, 2nd , 3rd, ... largest
% overlap. The assignment is bijective.
%
% INPUT:
% ======
% U_ref (required):
%   reference dictionary (num_dims x num_atoms)
%
% U (required):
%   dictionary (num_dims x num_atoms)
%
% OUTPUT:
% =======
% sim_vec:
%   vector containing the overlaps of matched atoms
%
% assign_vec:
%   defines how columns of U have to reordered to match witch column of
%   U_ref
function [sim_vec, assign_vec] = basis_similarity(U_ref, U)

num_dims = size(U, 1);
G = abs(U_ref'*U);

assign_vec = zeros(1,num_dims);
sim_vec = zeros(1,num_dims);

for cnt = 1:num_dims
    [~, idx] = max(G(:));
    [i,j] = ind2sub(size(G), idx);
    
    assign_vec(i) = j;
    sim_vec(cnt) = G(i, j);

    G(i,:) = 0;
    G(:,j) = 0;
end