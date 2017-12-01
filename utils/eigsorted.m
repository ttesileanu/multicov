function [V, D] = eigsorted(X, k, varargin)
% EIGSORTED Get eigenvalues and eigenvectors, sorted in decreasing order of
% eigenvalues.
%   E = EIGSORTED(X) returns the eigenvalues of the matrix in descending
%   order.
%
%   E = EIGSORTED(X, k) returns the top k eigenvalues of X.
%
%   [V, D] = EIGSORTED(X) returns the eigenvectors (the columns of V) and
%   the eigenvalues of the matrix X, in decreasing order of the eigenvalues.
%   The function also chooses the arbitrary sign of the eigenvectors such
%   that the largest element in absolute value is positive.
%
%   [V, D] = EIGSORTED(X, k) returns only the top k eigenvalues and
%   eigenvectors.

% XXX leaving the 'eigs' option out for now, since it doesn't really seem
% to help for non-sparse matrices

%   EIGSORTED(X, k, 'eigs') uses the Matlab eigs function instead of eig.
%   This can be much faster for large matrices if k is much smaller than
%   the size of the matrix.
%
% See also: EIG, EIGS.

% Tiberiu Tesileanu (2012, 2014)

% if no number is given, output all eigenvalues/eigenvectors
if nargin < 2
    k = size(X, 1);
end
% handle the 'eigs' option
if ~isempty(varargin)
    if ~strcmp(varargin{1}, 'eigs')
        error([mfilename ':badopt'], 'The third argument, if provided, should be ''eigs''.');
    end
    useeigs = true;
else
    useeigs = false;
end

if nargout < 2
    % only need eigenvalues
    if ~useeigs
        V = eig(X);
    else
        V = eigs(X, k);
    end
    V = sort(V, 'descend');
    V = V(1:k);
else
    % get eigenvalues and eigenvectors
    if ~useeigs
        [V, D] = eig(X);
    else
        [V, D] = eigs(X, k);
    end
    [D, indices] = sort(diag(D), 'descend');
    V = V(:, indices(1:k));
    D = D(1:k);

    for n = 1:k
        [~, idx] = max(abs(V(:, n)));
        % this can only be zero if V(:, n) is 0, which shouldn't be
        % possible... and if it is, multiplying it by 0 won't
        % change anything anyway
        factor = sign(V(idx(1), n));
        V(:, n) = factor*V(:, n);
    end
end

end