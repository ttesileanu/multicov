function [newalign, mapping] = alnresample(alignment, varargin)
% ALNRESAMPLE Resample sequences from alignment.
%   newalign = ALNRESAMPLE(alignment) generates a new alignment that is
%   obtained by uniform sampling of sequences, with replacement, from the
%   original alignment. The number of sequences in the new alignment is the
%   same as that in the old alignment. Because of the sampling with
%   replacement, some sequences will appear several times, and some will
%   not appear at all in newalign.
%
%   newalign = ALNRESAMPLE(alignment, n) makes an alignment of n sequences
%   instead of one that has the same size as the original alignment.
%
%   [newalign, mapping] = ALNRESAMPLE(...) also returns a vector describing
%   the mapping between the sequences in the new alignment and those in the
%   original one. 'Mapping' is a vector with as many entries as there are
%   sequences in newalign, and mapping(i) is the index in 'alignment' of
%   the sequence on the ith position in newalign.
%
%   Options:
%    'replacement' <b>
%       Set to false to perform sampling without replacement. In this case
%       of course n shouldn't be larger than the number of sequences in the
%       alignment.
%       (default: true)

% Tiberiu Tesileanu (2014-2015)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('n', size(alignment.data, 1), @(n) isnumeric(n) && isscalar(n));
parser.addParamValue('replacement', true, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

newalign = alignment;
if params.replacement
    mapping = randi([1 size(alignment.data, 1)], params.n, 1);
else
    mapping = randperm(size(alignment.data, 1), params.n);
end
newalign.data = alignment.data(mapping, :);
newalign.seqw = alignment.seqw(mapping);
newalign.annotations = alignment.annotations(mapping);

end