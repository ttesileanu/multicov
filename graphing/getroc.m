function [tp, fp] = getroc(prediction, truth, varargin)
% GETROC Calculate ROC curve for vector predictions.
%   [tp, fp] = GETROC(prediction, truth) calculates an ROC curve for the
%   given predictions, provided the true binary values.
%
%   Options:
%     'n' <n>
%       Maimum number of points to use. (default: 5000)

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('n', 5000, @(x) isscalar(x) && isnumeric(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

thresholds = sort(prediction(:));
truth = logical(truth);
n = min(length(thresholds), params.n);
tp = zeros(n+1, 1);
fp = zeros(n+1, 1);
% make sure we always start at 0, 0
for j = 1:n
    jj = floor(j*length(thresholds)/n);
    m = (prediction >= thresholds(jj));
    tp(j) = sum(m(:) & truth(:)) / sum(truth(:));
    fp(j) = sum(m(:) & ~truth(:)) / sum(~truth(:));
end

end