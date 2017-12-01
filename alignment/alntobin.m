function binalign = alntobin(alignment)
% ALNTOBIN Transform character alignment into binary alignment.
%   binalign = ALNTOBIN(alignment) creates a binary alignment containing
%   the information from the given character alignment. For each alphabet,
%   each position in a sequence is transformed into a list of 0s and 1s -- 1
%   at position A in this list means that the letter is the Ath one in the
%   alphabet. Gaps are not considered letters; they are coded by a list of
%   all 0s.
%
%   See also BINCHECK, BINTOALN.

% Tiberiu Tesileanu (2012-2014)

if ~alncheck(alignment)
    error([mfilename ':badarg'], 'First argument should be an alignment structure.');
end

% make sure the binary alignment has the right 'type' tag
binalign.type = 'binary';

nalpha = length(alignment.alphabets);
ranges = getalpharanges(alignment);

binalign.alphabets = alignment.alphabets;
binalign.alphawidths = alignment.alphawidths;
binalign.data = sparse([]);
for i=1:nalpha
    alphabet = alignment.alphabets{i};
    subdata = alignment.data(:, ranges(1, i):ranges(2, i));
    subbinalign = sparse(mattobin(subdata, alphabet));
    binalign.data = [binalign.data subbinalign];
end
binalign.seqw = alignment.seqw;
binalign.annotations = alignment.annotations;
binalign.refseq = alignment.refseq;

end