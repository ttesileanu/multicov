function res = crossvaldetails(fitfun, predfun, x, y, varargin)
% CROSSVALDETAILS Cross validation with detailed output.
%   res = CROSSVALDETAILS(FITFUN, PREDFUN, X, Y) performs cross-validation
%   when using FITFUN(XTRAIN, YTRAIN) to fit y as a function of x, assuming
%   that PREDFUN(PARAMS, XTEST) calculates the model predictions. PARAMS
%   here is the output from FITFUN. This calls Matlab's crossval function
%   internally.
%
%   Comapred to Matlab's crossval, CROSSVALDETAILS returns a structure that
%   contains both the outputs from the FUN function and predictions and
%   error estimates for the partitions.
%
%   All the options optional parameter name/value pairs from crossval are
%   supported. Additional options:
%    'evalfun':
%       A function that takes the model predictions and the actual expected
%       results, EVALFUN(PREDICTIONS, YTEST), and returns a scalar measure
%       of similarity between them.
%
%   See also: CROSSVAL.

n_params = [];

[local_opts, others] = splitoptions(varargin, {'evalfun'});

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('evalfun', [], @(f) isa(f, 'function_handle'));

% parse
parser.parse(local_opts{:});
opts = parser.Results;

    function out = fun(xtrain, ytrain, xtest, ytest)
        params = fitfun(xtrain, ytrain);
        predictions = predfun(params, xtest);
        errors = predictions - ytest;
        mean_error2 = mean(sum(errors.^2, 2), 1);
        
        if isempty(n_params)
            n_params = numel(params);
        elseif n_params ~= numel(params)
            error([mfilename ':badnparams'], 'Size of output from FITFUN should be the same for all partitions.');
        end
        
        if ~isempty(opts.evalfun)
            out = [params(:) ; mean_error2 ; opts.evalfun(predictions, ytest)];
        else
            out = [params(:) ; mean_error2];
        end
    end

cval_out = crossval(@fun, x, y, others{:});
res.params = cval_out(:, 1:n_params);
res.mse_each = cval_out(:, n_params+1);
res.mse = mean(res.mse_each);
if ~isempty(opts.evalfun)
    res.eval_each = cval_out(:, n_params+2);
    res.eval = mean(res.eval_each);
end

end