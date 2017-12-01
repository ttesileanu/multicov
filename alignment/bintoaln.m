function align = bintoaln(binalign)
% BINTOALN Transform binary alignment into character alignment.
%   align = BINTOALN(binalign) creates a character alignment containing
%   the information from the given binary alignment. For each alphabet,
%   each position in a sequence is transformed from a list of 0s and 1s to
%   a character -- 1 at position A in the binary list is transformed to the
%   Ath letter in the alphabet. Gaps are not considered letters; they are
%   coded by a list of all 0s.
%
%   See also ALNCHECK, ALNTOBIN.

% Tiberiu Tesileanu (2017)

if ~bincheck(binalign)
    error([mfilename ':badarg'], 'First argument should be a binary alignment structure.');
end

% make sure the binary alignment has the right 'type' tag
align.type = 'character';

nalpha = length(binalign.alphabets);
ranges = getalpharanges(binalign, 'binary');

align.alphabets = binalign.alphabets;
align.alphawidths = binalign.alphawidths;
align.data = '';
nseq = size(binalign.data, 1);
for i=1:nalpha
    alphabet = binalign.alphabets{i};
    letters = alphagetletters(alphabet);
    nlett = length(alphagetletters(alphabet, 'nogap'));
    subbindata = reshape(full(binalign.data(:, ranges(1, i):ranges(2, i)))', nlett, [], nseq);
    % to use the max below, need to do something to distinguish gaps from
    % the first alphabet letter...
    subbindata(nlett+1, :, :) = eps;
    
%    subidxdata = squeeze(sum(bsxfun(@times, (1:nlett)', subbindata), 1));
    [~, subidxdata] = max(subbindata, [], 1);
    subidxdata = squeeze(subidxdata);
    mask = (subidxdata <= nlett);
    subalign = repmat(letters(1), size(subidxdata));
    subalign(mask) = letters(1 + subidxdata(mask));
    
    align.data = [align.data subalign'];
end
align.seqw = binalign.seqw;
align.annotations = binalign.annotations;
align.refseq = binalign.refseq;

end