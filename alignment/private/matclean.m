function mat = matclean(mat, alphabet, replacement)
% MATCLEAN Internal function.
% 
% Clean the alignment data matrix by replacing all non-alphabet
% characters with some replacement character (by default: the gap).
%   newmat = matclean(mat, alphabet) replaces all unrecognized characters
%   in mat with the alphabet's gap (assumed to be the first character of
%   the alphabet). Note that for gapped alphabets (whose names start with
%   the letters 'gap'), the mock gap is ignored.
%
%   newmat = matclean(mat, alphabet, replacement) replaces all unrecognized
%   characters with replacement.
%
% See also: ALPHAGETLETTERS.

% Tiberiu Tesileanu (2012-2014)

if length(alphabet) > 3 && strcmp(alphabet(1:3), 'gap')
    alphabet = alphabet(4:end);
end

letters = alphagetletters(alphabet);
% decide on a replacement character
if nargin < 3
    replacement = letters(1);
end

% XXX this can fail in a variety of ways: if mat contains any \0
% characters, or if it contains any non-ASCII characters...
map = repmat(replacement, 1, 256);
map(letters) = letters;

% split the matrix in blocks and apply the map to each
% this way we avoid converting the whole matrix to double, which would take
% a lot of memory
blocksize = 128;
n = size(mat, 1);
for i = 1:blocksize:n
    idxs = i:min(i + blocksize - 1, n);
    mat(idxs, :) = map(mat(idxs, :));
end

end