function cdata = matautocorr(data, range, varargin)
% MATAUTOCORR Calculate autocorrelation function for a vector quantity.
%   cdata = MATAUTOCORR(data, range) calculates an autocorrelation function
%   treating each row of the data matrix as a sample. The function uses
%   the displacements from the 'range' array. Note that the elements in
%   'range' must not exceed size(data, 1)/2 - 1.
%
%   If range is empty or not provided, the whole possible range is used
%   (that is, 1:size(data, 1)/2-1).
%
%   The function assumes periodic boundary conditions, i.e., that the data
%   loops around at the end (data(end+1, :) == data(1, :)). For  displacements
%   much smaller than the maximum, this shouldn't make much of a difference.

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('slow', false, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

if nargin < 2
    range = [];
end

Nmax = size(data, 1)/2 - 1;
if isempty(range)
    range = 1:Nmax;
end

if any(range > Nmax) || any(range < 0)
    error([mfilename ':badrange'], 'Range values should be between 0 and size(data, 1)/2-1.');
end

range = range(:);
mu = mean(data, 1);
datactr = full(bsxfun(@minus, data, mu));
N = size(datactr, 1);
sigma2 = sum(datactr(:).^2) / (N-1);
% do not divide by 0!
if sigma2 < eps
    % no variation means perfect autocorrelation
    cdata = ones(size(range));
    return;
end

if params.slow
    cdata = zeros(size(range));
    for i = 1:length(range)
        k = range(i);
        dots1 = sum(datactr(1:end-k, :).*datactr(1+k:end, :), 2);
        dots2 = sum(datactr(end-k+1:end, :).*datactr(1:k, :), 2);
        cdata(i) = (sum(dots1) + sum(dots2))/((N - k)*sigma2);
    end
else
    al = fft(datactr);
    cdata0 = sum(ifft(abs(al).^2), 2);
    cdata = cdata0(range+1) ./ ((N-range).*sigma2);
end

end