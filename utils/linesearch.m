function [xmin, stats] = linesearch(fct, x0, x1, varargin)
% LINESEARCH Perform a line search for a minimum of the function.
%   xmin = LINESEARCH(fct, x0, x1) searches for a minimum of the 1d function
%   fct starting at x0, and moving in the direction in which the function
%   appears to be decreasing, based on the relation between fct(x0) and
%   fct(x1). This uses a golden section search. It is assumed that there is
%   only one minimum of the function.
%
%   xmin = LINSEARCH(fct, x0, x1, x2) assumes that fct(x0) >= fct(x1) and
%   fct(x2) >= fct(x1).
%
%   [xmin, stats] = LINESEARCH(...) returns some convergence statistics.
%   Stats is a structure with fields:
%       'result':   this can be 'converged' or 'failed'
%       'fmin':     minimum value of function that was found
%       'niter':    number of iterations
%       'neval':    number of function evaluations
%
%   Options:
%    'maxiter' <n>
%       Maximum number of iterations.
%       (default: 100)
%    'tol' <x>
%       Goal precision for fct(x).
%       (default: 1e-8)

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('x2', [], @(x) isscalar(x) && isnumeric(x));

parser.addParamValue('maxiter', 100, @(x) isscalar(x) && isnumeric(x));
parser.addParamValue('tol', 1e-8, @(x) isscalar(x) && isnumeric(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

x2 = params.x2;
phi = (1 + sqrt(5))/2;
if isempty(x2)
    % need to find an appropriate triplet of points
    
    % making sure f(x1) < f(x0)
    f0 = fct(x0);
    f1 = fct(x1);
    if f1 > f0
        tmp = x0;
        x0 = x1;
        x1 = tmp;
        
        tmp = f0;
        f0 = f1; %#ok<NASGU>
        f1 = tmp;
    end
    
    % a very simple strategy: go in the direction of decreasing fct until
    % we find the turning point
    x2 = x1 + phi*(x1 - x0);
    f2 = fct(x2);
    stats.neval = 3;
    while f1 > f2
        x0 = x1;
        x1 = x2;
        f0 = f1; %#ok<NASGU>
        f1 = f2;
        x2 = x1 + phi*(x1 - x0);
        f2 = fct(x2);
        stats.neval = stats.neval + 1;
    end
    
%     % this is translated from C from the Numerical Recipes book
%     % (and oh! do they write horrible C!)
%     
%     x2 = x1 + phi*(x1 - x0);
%     f2 = fct(x2);
%     stats.neval = 3;
%     ratio_limit = 100;
%     while f1 > f2
%         r = (x1 - x0)*(f1 - f2);
%         q = (x1 - x2)*(f1 - f0);
%         u = x1 - (q*(x1 - x2) - r*(x1 - x0)) / (2*(q - r));
%         ulim = x1 + ratio_limit*(x2 - x1);
%         if (x1 - u)*(u - x2) > 0
%             fu = fct(u);
%             stats.neval = stats.neval + 1;
%             if fu < f2
%                 x0 = x1;
%                 x1 = u;
%                 f1 = fu;
%                 break;
%             elseif fu > f1
%                 x2 = u;
%                 break;
%             else
%                 u = x2 + phi*(x2 - x1);
%                 fu = fct(u);
%                 stats.neval = stats.neval + 1;
%             end
%         elseif (x2 - u)*(u - ulim) > 0
%             fu = fct(u);
%             stats.neval = stats.neval + 1;
%             if fu < f2
%                 x1 = x2;
%                 x2 = u;
%                 u = x2 + phi*(x2 - x1);
%                 f1 = f2;
%                 f2 = fu;
%                 fu = fct(u);
%                 stats.neval = stats.neval + 1;
%             end
%         elseif (u - ulim)*(ulim - x2) >= 0
%             u = ulim;
%             fu = fct(u);
%             stats.neval = stats.neval + 1;
%         else
%             u = x2 + phi*(x2 - x1);
%             fu = fct(u);
%             stats.neval = stats.neval + 1;
%         end
%         x0 = x1;
%         x1 = x2;
%         x2 = u;
%         f0 = f1;
%         f1 = f2;
%         f2 = fu;
%     end
else
    f1 = fct(x1);
    stats.neval = 1;
end

% this is essentially a translation of Wikipedia's Java code
% with a different stopping criterion
%tau = sqrt(params.tol);
resphi = 2 - phi;
stats.result = 'failed';
lastfx = [];
% we're assuming it's true that fct(x0) > fct(x1) and fct(x2) > fct(x1)
for i = 1:params.maxiter
    if abs(x2 - x1) > abs(x1 - x0)
        x = x1 + resphi*(x2 - x1);
    else
        x = x1 + resphi*(x0 - x1);
    end
%     if abs(x2 - x0) < tau*(abs(x1) + abs(x))
%         stats.result = 'converged';
%         break;
%     end
    fx = fct(x);
    if ~isempty(lastfx) && abs(fx - lastfx) < params.tol
        stats.result = 'converged';
        break;
    end
    stats.neval = stats.neval + 1;
    if fx < f1
        if abs(x2 - x1) > abs(x1 - x0)
            x0 = x1;
        else
            x2 = x1;
        end
        x1 = x;
        f1 = fx;
    else
        if abs(x2 - x1) > abs(x1 - x0)
            x2 = x;
        else
            x0 = x;
        end
    end
    lastfx = fx;
end
stats.niter = i;
stats.fmin = f1;
xmin = x1;

end