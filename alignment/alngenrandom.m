function alignment = alngenrandom(N, freq, pattern)
% ALNGENRANDOM Generate random alignment.
%   alignment = ALNGENRANDOM(N, freq, pattern) generates an alignment of N
%   random sequences, drawn from the per-site probability distribution
%   identified by freq, and using the alphabet structure given by pattern.
%   Pattern must be a structure having fields alphabets and alphawidths,
%   identifying the alphabets used in the alignment. The freq vector has an
%   entry for each possible character (gap frequency should not be
%   contained in freq, except for gapped alphabets) for each position. For
%   example, for a protein alignment of length n, freq would have 20*n
%   elements.
%
%   If pattern contains a 'refseq' field, it is copied over to the
%   alignment.
%
%   ALNGENRANDOM(N, freq, alphabet) generates single-alphabet sequences
%   using the given alphabet. The size of the alignment is inferred from
%   the size of the frequency vector.
%
%   ALNGENRANDOM(N, pattern) generates sequences following the given
%   pattern using a uniform distribution for the characters within each
%   alphabet (gap included). In this form, the generated data is the same
%   whether alphabets are gapped or not.
%
%   ALNGENRANDOM(N, n, alphabet), where n is an integer, uses a uniform
%   distribution of characters to generate sequences of length n in the
%   given alphabet. In this form also, the generated data is the same
%   whether the alphabet is gapped or not.
%
%   See also: ALNMAKE, ALPHAGETLETTERS.

% Tiberiu Tesileanu (2012-2014)

if nargin == 2
    pattern = freq;
    if ~isstruct(pattern)
        error([mfilename ':badpattern'], 'In alngenrandom(N, pattern), the pattern should be a structure.');
    end
    freq = [];
end
if ischar(pattern)
    alphabet = pattern;
    clear pattern;
    pattern.alphabets = {alphabet};
    if isscalar(freq)
        n = freq;
        freq = [];
    else
        n = length(freq) / length(alphagetletters(alphabet, 'nogap'));
    end
    pattern.alphawidths = n;
end
if isempty(freq)
    for i = 1:length(pattern.alphabets)
        alphabet = pattern.alphabets{i};
        % for gapped alphabets, do not generate characters of the mock gap
        % variety
        if length(alphabet) > 3 && strcmp(alphabet(1:3), 'gap')
            alphabet_ng = alphabet(4:end);
        else
            alphabet_ng = alphabet;
        end
        nletts0 = length(alphagetletters(alphabet_ng));
        nletts = length(alphagetletters(alphabet, 'nogap'));
        freq = [freq ; ones(pattern.alphawidths(i)*nletts, 1)/nletts0]; %#ok<AGROW>
    end
end
width = sum(pattern.alphawidths);
data = repmat(blanks(width), N, 1);
nalpha = length(pattern.alphabets);
ranges = getalpharanges(pattern);
binranges = getalpharanges(pattern, 'binary');
for i = 1:nalpha
    alphabet = pattern.alphabets{i};
    
    subfreq = freq(binranges(1, i):binranges(2, i));
    data(:, ranges(1, i):ranges(2, i)) = matgenrandom(N, subfreq, alphabet);
end

% make sure to have all the right tags for an alignment
alignment = alnmake;

% fill out the members
alignment.alphabets = pattern.alphabets;
alignment.alphawidths = pattern.alphawidths;
alignment.data = data;
alignment.seqw = ones(N, 1);
alignment.annotations = repmat({''}, N, 1);
if ~isfield(pattern, 'refseq')
    for i = 1:nalpha
        alignment.refseq(i).seqdb = 'generic';
        alignment.refseq(i).seqid = '';
        alignment.refseq(i).map = (1:alignment.alphawidths(i))';
    end
else
    alignment.refseq = pattern.refseq;
end

end