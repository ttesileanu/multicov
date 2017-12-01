function cseq = getconsensus(alignment, varargin)
% GETCONSENSUS Calculate consensus sequence.
%   cseq = GETCONSENSUS(alignment) returns the consensus sequence for the
%   (character or binary) alignment. The consensus residue at position i is
%   the most common residue in this column. By default, this excludes gaps.
%   Note that it is irrelevant whether alphabets are gapped or not (i.e.,
%   whether their names start with 'gap').
%
%   GETCONSENSUS(..., 'gaps', true) takes gaps into consideration as
%   consensus residues.
%
%   GETCONSENSUS(stats, ...) uses a statistics structure instead of the
%   alignment. Note that only the freq1, alphabets, and alphawidths fields
%   are used.
%
%   Options:
%    'gaps' <b>
%       Whether to allow gaps to be considered consensus residues.
%       (default: false)
%    'indices' <b>
%       If this is true, it returns a vector with the indices of the
%       consensus residues in the binary alignment, as opposed to a
%       character array. Note that this does not work in conjunction with
%       'gaps', true.
%       (default: false)
%
% See also: GETFREQ1, GETSTATS, ALPHAGETLETTERS.

% Tiberiu Tesileanu (2012-2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('gaps', false, @(b) isscalar(b) && islogical(b));
parser.addParamValue('indices', false, @(b) isscalar(b) && islogical(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

if params.indices && params.gaps
    error([mfilename ':badcomb'], 'Incompatible combination of ''indices'', true and ''gaps'', true');
end

if statscheck(alignment)
    freq1 = alignment.freq1;
else
    freq1 = getfreq1(alignment);
end

ranges = getalpharanges(alignment, 'binary');
cseq = [];

for i = 1:length(alignment.alphabets);
    alphabet = alignment.alphabets{i};
    subfreq = freq1(ranges(1, i):ranges(2, i));
    if ~params.gaps
        letters = alphagetletters(alphabet, 'nogap');
        freq_block = reshape(subfreq, length(letters), alignment.alphawidths(i));
        if length(alphabet) > 3 && strcmp(alphabet(1:3), 'gap')
            freq_block = freq_block(2:end, :);
            letters = letters(2:end);
        end
    else
        letters = alphagetletters(alphabet);
        subfreq_reshaped = reshape(subfreq, length(letters) - 1, alignment.alphawidths(i));
        if length(alphabet) > 3 && strcmp(alphabet(1:3), 'gap')
            freq_block = subfreq_reshaped;
            letters = letters(2:end);
        else
            freq_block = [1 - sum(subfreq_reshaped, 1) ; subfreq_reshaped];
        end
    end
    [~, argmax] = max(freq_block, [], 1);
    if ~params.indices
        cseq = [cseq letters(argmax)]; %#ok<AGROW>
    else
        cseq = [cseq (ranges(1, i) + argmax + length(letters)*(0:(alignment.alphawidths(i)-1)) - 1)]; %#ok<AGROW>
    end
end

end