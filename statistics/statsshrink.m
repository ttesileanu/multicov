function newstats = statsshrink(stats, alpha, varargin)
% STATSSHRINK Shrink the statistics towards a background distribution.
%   newstats = STATSSHRINK(stats, alpha) shrinks the statistics towards a
%   background distribution. The amount parameter, alpha, is between 0 and 1, 0 meaning
%   no pseudocount (in which case stats is simply returned), 1 meaning none
%   of the original data is kept.
%
%   By default, this uses a uniform distribution of residues (including
%   gaps) for the background. See the 'background' option for alternatives.
%
%   Options:
%    'background' <s/v/c>
%       This changes the background distribution used by the function. Any
%       of the options allowed by extendfreq can be used.
%       (default: 'uniform')
%
% See also: GETDBBKG, EXTENDFREQ.

% Tiberiu Tesileanu (2013-2014)

if alpha < 0 || alpha > 1
    error([mfilename ':badalpha'], 'Alpha should be between 0 and 1.');
end

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

test_background = @(c) all(cellfun(@(s) ismember(s, {'uniform', 'database'}), c));
parser.addParamValue('background', 'uniform', ...
    @(c) (ischar(c) && test_background({c})) || (iscell(c) && test_background(c)));

% parse
parser.parse(varargin{:});
params = parser.Results;

newstats = stats;

% don't do anything if alpha is zero
if alpha > 0
    [bkgfreq1, bkgfreq2] = extendfreq(stats, params.background);
    bkgcmat = bkgfreq2 - bkgfreq1(:)*bkgfreq1(:)';
    
    newstats.cmat = (1 - alpha)*stats.cmat + alpha*bkgcmat;
    
    newstats.freq1 = (1 - alpha)*stats.freq1(:) + alpha*bkgfreq1(:);
    if isfield(stats, 'freq2')
        newstats.freq2 = newstats.cmat + newstats.freq1*newstats.freq1';
    end
end

end