% U = solve_orth_procrustes(M)
%
% Solves the orthogonal procrustes problem, a matrix approximation problem
% which asks for the orthogonal matrix U which is nearest to a given matrix
% M in terms of the Frobenius metric. If M = A'B, then the problem can be 
% equivalently formulated as follows:  U = arg min_W ||WA - B||_F subject 
% to W'W = I.
%
% [1] Schonemann, P.H. (1966), "A generalized solution of the orthogonal 
% Procrustes problem", Psychometrika 31: 1–10
%
% [2] https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem
%
% INPUT:
% ======
% M (required):
%   a matrix whose number of columns does not exceed it number of rows
%
% OUTPUT:
% =======
% U:
%   matrix with same size as M and mutually orthogonal columns.
%   Furthermore, U is closest to M in terms of the Frobenius metric

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
function U = solve_orth_procrustes(M)

[num_rows, num_columns] = size(M);
assert(num_columns <= num_rows);

[Q, S, R] = svd(M);

% set all singular values to one
S = eye(size(S));

U = Q*S*R';