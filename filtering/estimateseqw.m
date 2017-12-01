function seqw = estimateseqw(alignment, max_seqid, varargin)
% ESTIMATESEQW Estimate sequence weights for the alignment.
%   seqw = ESTIMATESEQW(alignment, max_seqid) estimates sequence weights by
%   treating sequences that are more than max_seqid similar as identical.
%   Similarity is calculated as the fraction of characters that match. More
%   precisely, the sequence weight for sequence A is 1 / n(A), where n(A)
%   is the number of sequences in the alignment that are more than max_seqid
%   identical to A. Preexisting sequence weights in the alignment are
%   ignored.

% Tiberiu Tesileanu (2012-2014)

if nargin == 1 && ischar(alignment) && strcmp(alignment, 'hascpp')
    seqw = (exist('calcseqw_cpp', 'file') == 3);
    return;
end

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('nocpp', false, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

%options = {'update', 'dist', 'no-cpp'};
%params = optparse(varargin, options);

if ~alncheck(alignment)
    error([mfilename ':badarg'], 'The first argument should an alignment structure.');
end

if ~params.nocpp && exist('calcseqw_cpp', 'file') == 3
    tmp = calcseqw_cpp(alignment.data, max_seqid);
    counts = 1 + sum(tmp + tmp');
else
    % since we're using Hamming distance, it doesn't really matter how we
    % convert the data to numbers
    numdata = double(alignment.data);
    counts = 1 + sum(squareform(pdist(numdata, 'hamm') < 1 - max_seqid), 1);
end

seqw = 1 ./ counts;

end