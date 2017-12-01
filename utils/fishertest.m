function p = fishertest(M, varargin)
% FISHERTEST Fisher Exact test for 2x2 contingency tables.
%   p = FISHERTEST(M) performs the non-parametric Fisher exact probability
%   test on a 2-by-2 contingency table described by M, and returns the
%   p-value. It calculates the exact probability of observing the given
%   and more extreme distributions on two variables.
%
%   Available options:
%     'tails' <n>
%       The argument can be 1 or -1, for the one-tailed test, or 2 for the
%       two-tailed one. See below for the difference between 1 and -1.
%       (default: 1)
%     'type' <s>
%       Options are:
%         'rel': the tail is chosen based on the data: if 'tails' is 1, the
%                tail is chosen in the direction of decreasing probability;
%                the opposite tail is chosen for 'tails' == -1.
%         'abs': the tail is chosen in the direction of increasing M(1, 1)
%                if 'tails' is 1, and in the direction of decreasing it if
%                'tails' is -1.
%       (default: 'rel')

% The function is adapted from Jos van der Geest's Mathworks contribution,
% http://www.mathworks.com/matlabcentral/fileexchange/23038-fishertest
%
% -------------------------------
% from Jon van der Geest's file:
% Source: Siegel & Castellan, 1988, "Nonparametric statistics for the
%         behavioral  sciences", McGraw-Hill, New York
%
% Created for Matlab R13+
% version 1.0 (feb 2009)
% (c) Jos van der Geest
% email: jos@jasen.nl
%
% File history:
%   1.0 (feb 2009) - created

% Adapted by Tiberiu Tesileanu (2013-2014)
%   The original code sometimes chose the wrong tail for one-tailed tests,
%   and lacked any support for two-tailed tests. The interface was also
%   changed to conform to MultiCOV's design.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('tails', 1, @(x) ismember(x, [1, -1, 2]));
parser.addParamValue('type', 'rel', @(s) ismember(s, {'rel', 'abs'}));

% parse
parser.parse(varargin{:});
params = parser.Results;

if ~isnumeric(M) || ~isequal(size(M), [2 2]) || any(M(:) ~= fix(M(:))) || any(M(:) < 0)
    error([mfilename ':badtable'], 'For numerical input, M should be a 2-by-2 matrix of positive integers.');
end

% The row and column sums are always the same
src = [sum(M, 2).' sum(M, 1)]; % row and column sums
N = src(1) + src(2);           % total number of obervations

% get the probability for the initial configuration
p0 = fratio(src, [M(:) ; N]);

dM = [1 -1 ; -1 1];

if abs(params.tails) == 1
    if strcmp(params.type, 'rel')
        % check which direction the tail points
        if M(1, 1)*M(2, 2) < M(1, 2)*M(2, 1)
            dM = -dM;
        end
    end
    
    if params.tails < 0
        dM = -dM;
    end
    
    steps = min(M(dM < 0));
    
    % now calculate the probality for observing the matrix M and all the
    % matrices that have more extreme distributions.
    p = p0;
    for i = 1:steps
        % make M more extreme
        M = M + dM;
        % calculate P value
        crtp = fratio(src, [M(:) ; N]);
        p = p + crtp;
    end
else
    % need tables on both tails
    % start with the most extreme table, and push to the other extreme
    M = M + dM*min(M(dM == -1));
    
    dM = -dM;
    p = 0;
    for i = 0:min(M(dM == -1))
        crtp = fratio(src, [M(:) ; N]);
        if crtp <= p0 + eps
            p = p + crtp;
        end
        M = M + dM;
    end
end

end

function fr = fratio(A, B)
% See FACTORIALRATIO for detailed help and comments
% http://www.mathworks.com/matlabcentral/fileexchange/23018
A = A(A>1) - 1;
B = B(B>1) - 1;
maxE = max([A(:) ; B(:) ; 0]);
if maxE > 1 && ~isequal(A, B)
    R = sparse(A, 1, 1, maxE, 2) + sparse(B, 2, 1, maxE, 2);
    R = flipud(cumsum(flipud(R)));
    R = (R(:, 1) - R(:, 2)).';
    X = 2:(maxE + 1);
    q = find(R);
    fr = prod(X(q) .^ R(q));
else
    fr = 1;
end

end