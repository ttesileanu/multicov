function ranges = getalpharanges(alignlike, opt)
% GETALPHARANGES Get the position ranges for each alphabet in an alignment.
%   ranges = GETALPHARANGES(alignment) returns a 2 x A matrix, where A is
%   the number of alphabets in the alignment, such that the ith alphabet
%   covers the columns in the range ranges(1, i):ranges(2, i).
%
%   Note that this only uses the 'alphawidths' field of an alignment, so
%   any structure with such a field can be used.
%
%   GETALPHARANGES(alphawidths) returns the result using just an array of
%   widths for each alphabet.
%
%   binranges = GETALPHARANGES(alignment, 'binary') returns the ranges in
%   binary alignment columns (each column of the usual alignment maps to L
%   columns in the binary alignment, where L is the number of characters in
%   the alignment *excluding* the gap). Note that this also uses the
%   'alphabets' field of the alignment.

% Tiberiu Tesileanu (2012-2014)

alphabets = {};
if isstruct(alignlike)
    % we're given an alignment-like structure
    alphawidths = alignlike.alphawidths;
    if isfield(alignlike, 'alphabets')
        alphabets = alignlike.alphabets;
    end
else
    alphawidths = alignlike;
end

if isempty(alphawidths) && isempty(alphabets)
    ranges = [];
    return;
end

if ~isnumeric(alphawidths) || ~isvector(alphawidths) || ~isreal(alphawidths)
    error([mfilename ':badwidth'], 'The argument should either be a numeric array or an alignment-like structure.');
end
if nargin >= 2 && strcmp(opt, 'binary')
    bin = true;
    if isempty(alphabets) && ~isempty(alphawidths)
        error([mfilename ':noalpha'], 'Need alphabet information to return binary alignment ranges.');
    end
else
    bin = false;
end

n = length(alphawidths);
ranges = zeros(2, n);
ranges(1, 1) = 1;
if ~bin
    ranges(2, n) = sum(alphawidths);
    for i = 1:n-1
        w = alphawidths(i);
        ranges(2, i) = ranges(1, i) + w - 1;
        ranges(1, i + 1) = ranges(2, i) + 1;
    end
else
    for i = 1:n
        w = alphawidths(i)*length(alphagetletters(alphabets{i}, 'nogap'));
        ranges(2, i) = ranges(1, i) + w - 1;
        if i < n
            ranges(1, i + 1) = ranges(2, i) + 1;
        end
    end
end

end