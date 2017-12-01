function [bool, details] = statscompare(stats1, stats2, varargin)
% STATSCOMPARE Compare two statistics structures.
%   bool = STATSCOMPARE(stats1, stats2) compares the two statistics
%   structures, and returns true if they are identical (within an error
%   tolerance; see options below).
%
%   [bool, details] = STATSCOMPARE(stats1, stats2) also returns a structure
%   detailing the comparison. This structure contains fields 'dalphabets'
%   and 'dalphawidths' that are true if the 'alphabets' or 'alphawidths'
%   fields are different between stats1 and stats2. If they are not
%   different, 'details' contains for each of the fields that are common
%   between the two structures ('freq1', 'cmat', and 'freq2' [if it
%   exists]), a field ('dfreq1', 'dcmat', 'or 'dfreq2') showing the
%   extent of the difference. Each of these fields is a two-element vector
%   in which the first element is the maximum absolute difference between
%   the fields of stats1 and stats2, and the second element is the mean
%   absolute difference.
%
%   If either of the statistics structures contains NaNs, the result is
%   false, even if the two structures are identical.
%
%   Options:
%    'tol' <x>
%       Allowed tolerance for the maximum deviation between fields of
%       stats1 and stats2.
%       (default: eps)

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('tol', eps, @(x) isnumeric(x) && isscalar(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

if ~statscheck(stats1) || ~statscheck(stats2)
    error([mfilename ':badstats'], 'The two input arguments should be statistics structures.');
end

details.dalphabets = ~isequal(stats1.alphabets(:), stats2.alphabets(:));
details.dalphawidths = ~isequal(stats1.alphawidths(:), stats2.alphawidths(:));
if any(isnan(stats1.freq1)) || any(isnan(stats2.freq1)) || ...
        any(isnan(stats1.cmat(:))) || any(isnan(stats2.cmat(:))) || ...
        (isfield(stats1, 'freq2') && any(isnan(stats1.freq2(:)))) || ...
        (isfield(stats2, 'freq2') && any(isnan(stats2.freq2(:))))
    bool = false;
    details.dfreq1 = [nan nan];
    details.dcmat = [nan nan];
elseif ~details.dalphabets && ~details.dalphawidths && ...
        isequal(size(stats1.freq1), size(stats2.freq1)) && ...
        isequal(size(stats1.cmat), size(stats2.cmat))
    difffreq1 = abs(stats1.freq1 - stats2.freq1);
    diffcmat = abs(stats1.cmat - stats2.cmat);
    details.dfreq1 = [max(difffreq1) mean(difffreq1)];
    details.dcmat = [max(diffcmat(:)) mean(diffcmat(:))];
    if details.dfreq1(1) > params.tol || details.dcmat(1) > params.tol
        bool = false;
    else
        bool = true;
    end
    if isfield(stats1, 'freq2') && isfield(stats2, 'freq2')
        difffreq2 = abs(stats1.freq2 - stats2.freq2);
        details.dfreq2 = [max(difffreq2(:)) mean(difffreq2(:))];
        if details.dfreq2(1) > params.tol
            bool = false;
        end
    end
    if isfield(stats1, 'freq2') ~= isfield(stats2, 'freq2')
        bool = false;
    end
else
    bool = false;
end

end