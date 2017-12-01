function [p, chi2] = chi2gof2(m, varargin)
% CHI2GOF2 Perform a multi-sample chi-squared test.
%   p = CHI2GOF2(m) returns the p-value for a multi-sample chi-squared
%   test, where the n x k matrix m is the contingency table. The null
%   hypothesis can be interpreted as the statement that each of the k
%   samples represented by the k columns of m are drawn from the same
%   population. Note, however, that the test is invariant under transposing
%   the matrix m, so the interpretation is the same if columns and rows are
%   exchanged. This is a two-tailed test.
%
%   [p, chi2] = CHI2GOF2(m) also returns the value of the chi^2 statistic.
%
%   Options:
%    'yates' <b>
%       Whether to use the Yates correction for continuity.
%       (default: false in general, true if m is 2 x 2)

% Tiberiu Tesileanu (2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('yates', [], @(b) isscalar(b) && islogical(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

if isempty(params.yates)
    params.yates =  all(size(m) == [2 2]);
end

sums1 = sum(m, 1);
sums2 = sum(m, 2);
expected = sums2*sums1/sum(sums1);
diff = m - expected;

if params.yates
    diff = abs(diff) - 0.5;
end

% make sure to avoid dividing by zero
ratio = zeros(size(diff));
mask = (expected ~= 0);
ratio(mask) = diff(mask).^2 ./ expected(mask);
chi2 = sum(ratio(:));

df = prod(size(m) - 1);
p = 1 - chi2cdf(chi2, df);

end