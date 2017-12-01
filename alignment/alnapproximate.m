function [newalign, maps] = alnapproximate(alignment, n, varargin)
% ALNAPPROXIMATE Make an approximation of the alignment by keeping only the
% top most common characters in each column.
%   newalign = ALNAPPROXIMATE(alignment, n) creates an n-letter
%   approximation of the alignment. This is done by keeping only the top
%   n-1 most common characters in each column, and replacing everything else
%   with '0'. The top n-1 characters are encoded by '1', '2', etc., in the
%   order of their abundance (the caracters used are the ones employed by
%   the appropriate 'multi' alphabet; see alphagetletters). Note that the
%   mapping from the 'multi' alphabets to the original ones usually differs
%   from column to column. An error is generated if any of the alphabets
%   used in the alignment contains less than n characters (including the
%   gap).
%
%   If the alignment has non-trivial sequence weights, these are taken into
%   account when estimating the abundance of characters. The result
%   essentially coincides with the binary approximation when n = 2.
%
%   [newalign, maps] = ALNAPPROXIMATE(...) also returns a character
%   matrix, maps, with n-1 rows. Each column contains the top n-1 characters
%   found in the corresponding alignment column, in decreasing order of
%   abundance. This matrix can be used to partially map back from the
%   approximate alignment to the original one.
%
%   Here are some options that can be provided to the function:
%    'gapszero' <b>
%       Whether to treat gaps separately by not allowing them to be
%       included in the top n characters. When this option is provided,
%       gaps are always mapped to '0'.
%       (default: true)
%
% See also: ALPHAGETLETTERS, ALNTOBIN.

% Tiberiu Tesileanu (2013-2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('gapszero', true, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

if n < 2 || n > 36
    error([mfilename ':badn'], 'Need 2 <= n <= 36.');
end

newalign = alignment;

nalphas = length(alignment.alphabets);
Nseqs = size(alignment.data, 1);
ranges = getalpharanges(alignment);

multiname = ['multi' int2str(n)];
multiletters = alphagetletters(multiname);

if nargout > 1
    maps = repmat(blanks(size(alignment.data, 2)), n-1, 1);
end
% XXX do we want to preserve the old alphabet information?
% newalign.consalphabets = newalign.alphabets;
% newalign.consalphawidths = newalign.alphawidths;
% newalign.consensus = blanks(size(alignment.data, 2));
for a = 1:nalphas
    alphabet = alignment.alphabets{a};
    allletters = alphagetletters(alphabet);
    if params.gapszero
        letters = alphagetletters(alphabet, 'nogap');
    else
        letters = allletters;
    end
    if length(allletters) >= n
        newalign.alphabets{a} = multiname;
        for i = ranges(1, a):ranges(2, a)
            % find the top n-1 most common characters in this column
            column = alignment.data(:, i);
            counts = zeros(1, length(letters));
            for j = 1:length(letters)
                lett = letters(j);
                counts(j) = alignment.seqw(:)'*(column == lett);
            end
            [~, order] = sort(counts, 'descend');
            tokeep = letters(order(1:n-1));
            if nargout > 1
                maps(:, i) = tokeep;
            end
            
%             newalign.consensus(i) = tokeep(1);
            
            % replace by the 'multi' alphabet
            for j = 1:n-1
                % start by making everything equal to the 'other'
                % character, the one that represents stuff that's not in
                % the top n
                newalign.data(:, i) = repmat(multiletters(1), Nseqs, 1);
                % then do the mapping for the top n characters
                for k = 1:n-1
                    mask = (column == tokeep(k));
                    newalign.data(mask, i) = multiletters(k + 1);
                end
            end
        end
    else
        error([mfilename ':smallalpha'], 'All alphabets must have at least n characters (including the gap).');
    end
end

end
