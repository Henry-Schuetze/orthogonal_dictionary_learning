% [sim_vec, match_atoms_vec] = dictionary_similarity(C, C_ref)
%
% Computes the overlaps between "matching columns" of matrix C and matrix 
% C_ref. The matching is assigned as a bijective mapping based on 1st, 2nd,
% 3rd, ... largest overlap. The columns of C and C_ref should all have unit
% length.
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
%   similarity vector (num_atoms x 1) containing the overlaps of matched
%   atom pairs
%
% match_atoms_vec:
%   row vector that matches the columns of matrix C to those of C_ref.
%   defines column permutation of matrix C such that a column
%   C(:, match_atoms_vec(col_idx)) matches with column C_ref(:, col_idx)
%

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
function [sim_vec, match_atoms_vec] = dictionary_similarity(C, C_ref)

% assert(isequal(size(C), size(C_ref)), 'C and C_ref do not have equal size')

num_atoms = size(C, 2);
G = abs(C_ref'*C);

match_atoms_vec = zeros(1, num_atoms);
sim_vec = zeros(num_atoms, 1);

for h = 1:num_atoms
    [~, idx] = max(G(:));
    [i,j] = ind2sub(size(G), idx);
    
    match_atoms_vec(i) = j;
    sim_vec(h) = G(i, j);

    G(i,:) = 0;
    G(:,j) = 0;
end