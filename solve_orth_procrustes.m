% U = solve_orth_procrustes(M)
%
% Solves the orthogonal procrustes problem, a matrix approximation problem
% which asks for the orthogonal matrix U which is nearest to a given matrix
% M in terms of Frobenius metric. Assuming M = A'B, it can be equivalently 
% formulated as follows:  U = arg min_W ||WA - B||_F subject to W'W = I.
%
% [1] Schonemann, P.H. (1966), "A generalized solution of the orthogonal 
% Procrustes problem", Psychometrika 31: 1–10
%
% [2] https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem
%
% INPUT:
% ======
% M (required):
%   a given matrix
%
% OUTPUT:
% =======
% U:
%   matrix of same size as M with orthogonal columns which is closest to M
%   in terms of Frobenius metric

% Henry Schuetze 
% Institute for Neuro- and Bioinformatics
% University of Luebeck, Germany
% Henry.Schuetze@uni-luebeck.de
function U = solve_orth_procrustes(M)

[num_dims, num_atoms] = size(M);
assert(num_atoms <= num_dims);

[Q, S, R] = svd(M);

% set all singular values to one
S = eye(size(S));

U = Q*S*R';