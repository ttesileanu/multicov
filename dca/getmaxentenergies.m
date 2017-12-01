function energies = getmaxentenergies(alignment, maxentparams, varargin)
% GETMAXENTENERGIES Calculate the maximum entropy model energies associated
% with a set of sequences.
%   energies = GETMAXENTENERGIES(alignment, maxentparams) uses the maximum
%   entropy model parameters to calculate energies for each of the
%   sequences in the alignment. The structure of the alignment must match
%   that of maxentparams (though 'gapped' and 'non-gapped' alphabets are
%   treated as the same). The alignment can be a character or binary
%   alignment.
%
%   energies = GETMAXENTENERGIES(seqmat, maxentparams) calculates energies
%   for the sequences given by the caracter matrix seqmat. The number of
%   columns of this matrix should match the alphawidths from maxentparams.

% Tiberiu Tesileanu (2014)

hascpp = (exist('cppevalmaxent', 'file') == 3);
if nargin == 1 && ischar(alignment) && strcmp(alignment, 'hascpp')
    energies = hascpp;
    return;
end

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('nocpp', false, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

maxentparams = paramsincludegaps(maxentparams);
if ~params.nocpp && hascpp
    isaln = alncheck(alignment);
    if isaln || (ischar(alignment) && ismatrix(alignment))
        if isaln
            data = alignment.data;
        else
            data = alignment;
        end
        alphaletters = cellfun(@(a) alphagetletters(a, 'nogap'), maxentparams.alphabets, 'uniform', false);
        energies = cppevalmaxent(data', alphaletters, double(maxentparams.alphawidths), double(maxentparams.couplings));
        return;
    end
end
if alncheck(alignment)
    binalign = alntobin(alnincludegaps(alignment));
elseif bincheck(alignment)
    binalign = alnincludegaps(alignment);
elseif ischar(alignment) && ismatrix(alignment)
    alignment = alnincludegaps(struct(...
        'type', 'character', 'alphabets', {maxentparams.alphabets}, ...
        'alphawidths', {maxentparams.alphawidths}, ...
        'refseq', {maxentparams.refseq}, ...
        'annotations', {{}}, ...
        'data', {alignment}, 'seqw', []));
    binalign = alntobin(alignment);
else
    error([mfilename ':badarg'], 'The first argument should be an alignment or a character matrix.');
end

energies = -(1/2)*full(dot(binalign.data', maxentparams.couplings*binalign.data'));

end