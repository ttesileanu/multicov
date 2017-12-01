function newalign = filterrows(alignment, varargin)
% FILTERROWS Filter the alignment, eliminating sequences with too many
% gaps or sequences that are too far from a reference sequence.
%   newalign = FILTERROWS(alignment, 'gaps', max_gaps) eliminates sequences
%   that have a fraction of gaps larger than max_gaps.
%
%   newalign = FILTERROWS(alignment, 'minrefsim', min_seqid) eliminates
%   sequences that have less than min_seqid identity to the reference
%   sequence, which is defined here as the first sequence in the alignment.
%
%   If both 'gaps' and 'minrefsim' are given, the gapped sequences are
%   eliminated first (regardless of the order in which the flags are used).

% Tiberiu Tesileanu (2012-2013)

if ~alncheck(alignment)
    error([mfilename ':badarg'], 'The first argument should an alignment structure.');
end

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('gaps', 1, @(x) isnumeric(x) && isscalar(x));
parser.addParamValue('minrefsim', 0, @(x) isnumeric(x) && isscalar(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

if params.gaps < 1
    gaps = getgapstructure(alignment);
    
    % get rid of sequences with too many gaps
    gapfraction = sum(gaps, 2) / size(gaps, 2);
    mask1 = (gapfraction <= params.gaps);
else
    mask1 = true(1, size(alignment.data, 1));
end

newdata = alignment.data(mask1, :);
newseqw = alignment.seqw(mask1);
newannot = alignment.annotations(mask1);

if params.minrefsim > 0
    % calculate similarities to reference sequence
    refseq = alignment.data(1, :);
    simmat = (newdata == repmat(refseq, size(newdata, 1), 1));
    sim = sum(simmat, 2) / size(simmat, 2);
    
    % get rid of sequence too far from reference
    mask2 = (sim >= params.minrefsim);
else
    mask2 = true(1, size(newdata, 1));
end

% make sure to copy all the fields of the input structure
newalign = alignment;

newalign.data = newdata(mask2, :);
newalign.seqw = newseqw(mask2);
newalign.annotations = newannot(mask2);

end
