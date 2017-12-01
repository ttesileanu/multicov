function [me_params, minstats] = getexactmaxent(stats, varargin)
% GETEXACTMAXENT Find maximum entropy model that best matches the given
% statistics, using an iterative method.
%   me_params = GETEXACTMAXENT(stats) finds the maximum entropy model best
%   matching the given statistics using a slow but exact iterative method.
%
%   [me_params, minstats] = GETEXACTMAXENT(...) also returns some statistics
%   from the minimization procedure.
%
%   Options:
%    'dispinterval' <x>
%       How often to display progress information, in seconds.
%       (default: 10)
%    'function' <f>
%       Function to use for getting statistics from maximum entropy
%       parameters. This should be a handle for a function that will be
%       given the parameters structure and should return a statistics
%       structure.
%       (default: @getmaxentstatsslow without verbosity)
%    'maxsteps' <n>
%       Maximum number of steps to take.
%       (default: 10000)
%    'solver' <x>
%       Which solving method to use (see findminimum's 'method').
%       (default: 'cg')
%    'tol' <x>
%       Tolerance to use when comparing the target structure with the
%       current one.
%       (default: 1e-6)
%    'verbose' <b>
%       Whether to output progress information.
%       (default: true)

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('dispinterval', 10, @(x) isscalar(x) && isnumeric(x));
parser.addParamValue('function', @(p) getmaxentstatsslow(p, 'verbose', false), ...
    @(f) isscalar(f) && isa(f, 'function_handle'));
parser.addParamValue('maxsteps', 10000, @(n) isscalar(n) && isnumeric(n));
parser.addParamValue('tol', 1e-6, @(x) isscalar(x) && isnumeric(x));
parser.addParamValue('verbose', true, @(b) isscalar(b) && islogical(b));
parser.addParamValue('solver', 'cg', @(s) ismember(s, {'descent', 'cg'}));

% parse
parser.parse(varargin{:});
params = parser.Results;

% make sure we don't have the gaps in the statistics structure -- this is
% used to fix the gauge (we're fixing it to gap gauge)
stats = statsexcludegaps(stats);
% make sure we have the 'freq2' field
stats = getstats(stats, 'freq2', true);

% first we need a guess for me_params
me_params.type = 'maxent';
me_params.alphabets = stats.alphabets;
me_params.alphawidths = stats.alphawidths;
me_params.refseq = stats.refseq;

binmap = getbinmap(stats);
me_params.couplings = diag(2*cell2mat(...
    cellfun(@(idxs) flatten(log(stats.freq1(idxs) / (1 - sum(stats.freq1(idxs))))), ...
    binmap(:), 'uniform', false)));

% define the function to minimize, and its derivatives
fct = @(x) negloglikelihood(x, me_params, stats, params.function);
fctder = @(x) negloglikelihoodderivative(x, me_params, stats, params.function);

[x0, minstats] = findminimum(fct, triu(me_params.couplings), 'derivative', fctder, ...
    'verbose', params.verbose, 'dispinterval', params.dispinterval, ...
    'maxiter', params.maxsteps, 'tol', params.tol, 'method', params.solver);

me_params = paramsfromtriu(x0, me_params);

end

function [params, couplingstriu1] = paramsfromtriu(couplingstriu, structure)
% PARAMSFROMTRIU Recover parameters structure from triangular form.
%   [params, couplingstriu] = PARAMSFROMTRIU(couplingstriu, structure)
%   recovers a parameters structure given the upper triangular form
%   couplingstriu and a template structure.
%
%   The function also returns a version of couplingstriu with the diagonal
%   removed (couplingstriu1).

params = structure;
couplingstriu1 = couplingstriu - diag(diag(couplingstriu));
params.couplings = couplingstriu + couplingstriu1';

end

function nlL = negloglikelihood(couplingstriu, structure, stats0, fct)
% NEGLOGLIKELIHOOD The function to minimize, the negative log likelihood.
%   nlL = NEGLOGLIKELIHOOD(couplingstriu, structure) calculates the negative
%   log likelihood given the upper triangular portion of the coupling
%   matrix (couplingstriu), a parameters structure with the right
%   alphabets and alphabet sizes (structure), and a set of 'target'
%   statistics ('stats0'). The calculation is done using the function handle
%   fct.

[params, couplingstriu1] = paramsfromtriu(couplingstriu, structure);
[~, Z] = fct(params);
nlL = log(Z) - 0.5*dot(diag(couplingstriu), stats0.freq1) - ...
    dot(flatten(couplingstriu1), stats0.freq2(:));

end

function der = negloglikelihoodderivative(couplingstriu, structure, stats0, fct)
% NEGLOGLIKELIHOODDERIVATIVE The derivative of the objective function.
%   der = NEGLOGLIKELIHOODDERIVATIVE(couplingstriu, structure, stats0)
%   calculates the derivatives of the objective function with respect to the
%   parameters. The function is given the upper triangular portion of the
%   coupling matrix (couplingstriu1), a parameters structure with the right
%   alphabets and alphabet sizes (structure), and a set of 'target'
%   statistics 'stats0' (which should include the 'freq2' field). The
%   calculation is done using the function handle fct.

params = paramsfromtriu(couplingstriu, structure);
stats = fct(params);
der = 0.5*diag(stats.freq1 - stats0.freq1) + triu(stats.freq2 - stats0.freq2);

end