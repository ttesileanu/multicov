function newstats = statspseudocount(stats, alpha, varargin)
% STATSPSEUDOCOUNT Add pseudocount to the statistics.
%   newstats = STATSPSEUDOCOUNT(stats, alpha) adds pseudocount to the
%   statistics. The amount parameter, alpha, is between 0 and 1, 0 meaning
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
%    'fieldalpha' <x>
%       If provided, the function uses a different pseudocount amount for
%       the fields compared to the two-point function and correlations.
%       Note that this will break the relation cmat = freq2 - freq1*freq1'
%       for the resulting structure, newstats.
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
parser.addParamValue('fieldalpha', [], @(x) isnumeric(x) && isscalar(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

newstats = stats;

% don't do anything if alpha is zero
if alpha > 0
    [bkgfreq1, bkgfreq2] = extendfreq(stats, params.background);

    % it's easiest to add the pseudocount to freq2 instead of cmat
    if isfield(stats, 'freq2')
        freq2 = stats.freq2;
    else
        freq2 = stats.cmat + stats.freq1(:)*stats.freq1(:)';
    end
    freq2 = (1 - alpha)*freq2 + alpha*bkgfreq2;
    if isfield(stats, 'freq2')
        newstats.freq2 = freq2;
    end
    freq1_0 = (1 - alpha)*stats.freq1(:) + alpha*bkgfreq1(:);
    newstats.cmat = freq2 - freq1_0*freq1_0';
    
    if ~isempty(params.fieldalpha)
        newstats.freq1 = (1 - params.fieldalpha)*stats.freq1(:) + ...
            params.fieldalpha*bkgfreq1(:);
    else
        newstats.freq1 = freq1_0;
    end
end

end