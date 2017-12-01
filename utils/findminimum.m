function [xmin, stats] = findminimum(fct, x0, varargin)
% FINDMINIMUM Find minimum of (potentially non-quadratic) function.
%   xmin = FINDMINIMUM(fct, x0) searches for a minimum of the given
%   function, using x0 as a guess value. Currently the only search
%   algorithms implemented are gradient based, so the 'derivative' option
%   is mandatory. The dimensionality of the problem is given by the size
%   (and shape) of the guess value, x0.
%
%   Matlab's Optimization Toolbox has more powerful functions performing
%   minimization. Those should be used whenever possible. The purpose of
%   this function is to provide a reasonable replacement when the
%   Optimization Toolbox is not available.
%
%   [xmin, stats] = FINDMINIMUM(...) returns some statistics regarding the
%   search. It is a structure with fields
%       'result':   this can be 'converged' or 'failed'
%       'fmin':     minimum value of function that was found
%       'niter':    number of iterations
%       'nfeval':   number of function evaluations
%       'ndeval':   number of derivative evaluations
%
%   Options:
%    'derivative' <handle>
%       This is a function handle returning the derivative of the objective
%       function with respect to each of the variables (the gradient). The
%       result thus should have the same size as x0. Since all the algorithms
%       currently implemented are gradient-based, it is mandatory to
%       provide the derivative.
%    'dispinterval' <x>
%       How often to display progress information (in seconds).
%       (default: 10)
%    'maxiter' <n>
%       Maximum number of iterations to perform.
%       (default: 1e4)
%    'method' <s>
%       Algorithm to use:
%        'cg':          conjugate-gradient
%        'descent':     gradient descent
%       (default: 'cg')
%    'rate' <x>
%       Initial step to guess for line search (in units of the gradient).
%       (default: 1)
%    'tol' <x>
%       Value for the norm of the derivative below which convergence is
%       assumed.
%       (default: 1e-6)
%    'verbose' <b>
%       Whether to display progress information.
%       (default: true)

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('derivative', [], @(f) isa(f, 'function_handle'));
parser.addParamValue('dispinterval', 10, @(x) isscalar(x) && isnumeric(x));
parser.addParamValue('maxiter', 1e4, @(x) isscalar(x) && isnumeric(x));
parser.addParamValue('method', 'cg', @(s) ismember(s, {'descent', 'cg'}));
parser.addParamValue('rate', 1, @(x) isscalar(x) && isnumeric(x));
parser.addParamValue('tol', 1e-6, @(x) isscalar(x) && isnumeric(x));
parser.addParamValue('verbose', true, @(x) isscalar(x) && islogical(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

x = x0;
xmin = x;
stats.result = 'failed';
stats.fmin = fct(x);
stats.nfeval = 1;
stats.ndeval = 0;
prevder = [];
prevdirection = [];
timer = tic;
for i = 1:params.maxiter
    if params.verbose
        dt = toc(timer);
        if dt >= params.dispinterval
            disp(['findminimum iteration ' int2str(i) ', fmin ' num2str(stats.fmin) '...']);
            timer = tic;
        end
    end
    der = params.derivative(x);
    stats.ndeval = stats.ndeval + 1;
    switch params.method
        case 'descent'
            direction = der;
        case 'cg'
            if isempty(prevder)
                % first step is just gradient descent
                direction = der;
            else
                % Polak-Ribiere method
                beta = max(dot(der(:), der(:) - prevder(:)) / dot(prevder(:), prevder(:)), 0);
                direction = der + beta*prevdirection;
            end
    end
%     if norm(der(:))/numel(der) < params.tol
%         stats.result = 'converged';
%         break;
%     end
    [alpha, substats] = linesearch(@(a) fct(x + a*direction), 0, -params.rate, 'tol', params.tol);
    x = x + alpha*direction;
    stats.nfeval = stats.nfeval + substats.neval;
    f = substats.fmin;
    lastf = stats.fmin;
    if f < stats.fmin
        stats.fmin = f;
        xmin = x;
    end
    if abs(stats.fmin - lastf) < params.tol
        stats.result = 'converged';
        break;
    end
    if ~strcmp(substats.result, 'converged')
        break;
    end
    prevder = der;
    prevdirection = direction;
end
stats.niter = i;

end