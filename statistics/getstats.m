function stats = getstats(alignment, varargin)
% GETSTATS Calculate statistics of the alignment.
%   stats = GETSTATS(alignment) calculates positional frequencies and joint
%   parwise frequencies for the columns of the given alignment. This can be
%   either a character or a binary alignment.
%
%   The function returns a structure with the following fields
%    'freq1' -- the frequencies of residues at each site
%    'freq2' -- the frequencies of pairs of residues (this field is only
%               present if the 'freq2', true option is used)
%    'cmat'  -- the covariance matrix for the alignment
%    'alphabets', 'alphawidths'
%            -- copied from the alignment
%    'nseqs' -- number of sequences in the alignment
%    'type'  -- identifier of the type of structure; set to 'stats'
%
%   stats = GETSTATS(stats, ...) can be used to either add or remove a
%   'freq2' field from a statistics structure.
%
%   Options:
%    'freq2' <b>
%       Whether to also include the pairwise frequencies, or only the
%       covariance matrix. (default: false)
%
% See also: GETFREQ1, GETFREQ2, ALNTOBIN.

% Tiberiu Tesileanu (2013-2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('freq2', false, @(b) isscalar(b) && islogical(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

if alncheck(alignment)
    binalign = alntobin(alignment);
elseif bincheck(alignment)
    binalign = alignment;
elseif statscheck(alignment)
    stats = alignment;
    if params.freq2 && ~isfield(stats, 'freq2')
        stats.freq2 = stats.cmat + stats.freq1(:)*stats.freq1(:)';
    elseif ~params.freq2 && isfield(stats, 'freq2')
        stats = rmfield(stats, 'freq2');
    end
    return;
else
    error([mfilename ':badaln'], 'The argument must be an alignment structure, either binary or character, or a statistics structure.');
end

stats.type = 'stats';
stats.freq1 = getfreq1(binalign);
freq2 = getfreq2(binalign);
if params.freq2
    stats.freq2 = freq2;
end
stats.cmat = freq2 - stats.freq1(:)*stats.freq1(:)';
stats.alphabets = binalign.alphabets;
stats.alphawidths = binalign.alphawidths;
stats.refseq = binalign.refseq;
stats.annotations = binalign.annotations;
stats.nseqs = size(alignment.data, 1);

% XXX other fields?
%stats = copyfields(stats, binalign, {'consensus', 'consalphabets', 'consalphawidths', 'refseq'});

end