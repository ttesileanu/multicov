function reduced = reducecouplings(maxentparams, method, varargin)
% REDUCECOUPLINGS Reduce a maximum entropy model coupling matrix to a
% matrix of pairwise interaction strengths between residues.
%   reduced = REDUCECOUPLINGS(maxentparams, method) takes in the results
%   from a maximum entropy model estimation and calculates a matrix in
%   which the (i, j) entry gives an estimate for the strength of
%   interaction between residues i and j. The reduction is performed using
%   one of the following methods:
%    'di':   use the 'direct information' approach from Weigt et al. 2009
%            and Morcos et al. 2011; in this approach, reduced(i, j) is a
%            mutual information calculated for sites i and j using the
%            "direct" probability model estimated by calculate directmodel.
%            Note that this model requires the frequency vector for the
%            alignment for which the couplings were calculated; this needs
%            to be provided using the 'freq1' option.
%    'norm': take a norm of each sub-block of the matrix, as indicated by
%            the 'norm' option (see below); note that this is not invariant
%            to gauge changes. See the 'gauge' option for this.
%
%   Note that this function works on the version of the maximum entropy
%   parameters that includes gaps. If the input does not, paramsincludegaps
%   is called first.
%
%   Options:
%    'freq1' <v>
%       A vector of frequencies that needs to be provided for the 'di'
%       method.
%    'gauge' <s>
%       This is used to set the gauge in which the norm is calculated when
%       the method is 'norm'. It is ignored for 'di'. This can be
%        'original':    keep the gauge in which the parameters were provided
%        'zerosum':     perform the calculation in the zero-sum gauge
%       (default: 'zerosum')
%    'norm' <f/s>
%       Function handle for the norm function to be used when the method is
%       'norm'. This can also be a string, in which case it can be
%        'frobenius':   use the frobenius norm, @(m) norm(m, 'fro')
%        'spectral':    use the spectral norm (largest singular value), @norm
%       (default: 'frobenius')
%
% See also: CALCULATEDIRECTMODEL, CHANGEGAUGE, PARAMSINCLUDEGAPS.

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('freq1', [], @(v) isvector(v) && isnumeric(v));
parser.addParamValue('gauge', 'zerosum', @(s) isvector(s) && ischar(s));
parser.addParamValue('norm', 'frobenius', ...
    @(s) (ischar(s) && ismember(s, {'frobenius', 'spectral'})) || isa(s, 'function_handle'));

% parse
parser.parse(varargin{:});
params = parser.Results;

switch method
    case 'di'
        if isempty(params.freq1)
            error([mfilename ':missfreq1'], 'The ''freq1'' parameter needs to be set when the method is ''di''.');
        end
        dmodel = calculatedirectmodel(maxentparams, params.freq1);
        binmap = getbinmap(dmodel);
        % the function below looks like a mess, but it's faster to write
        % this in one function rather than define and call auxiliary
        % functions...
        reduced = blockapply(dmodel.freq2, @(m, i, j) ...
            trace(m' * log(m ./ ...
                (flatten(dmodel.freq1(binmap{i})) * ...
                 transpose(flatten(dmodel.freq1(binmap{j})))) ...
            )), dmodel, 'indices', true);
    case 'norm'
        % make sure our parameters include the gaps
        maxentparams = paramsincludegaps(maxentparams);
        if ischar(params.norm)
            switch params.norm
                case 'frobenius'
                    params.norm = @(m) norm(m, 'fro');
                case 'spectral'
                    params.norm = @norm;
                otherwise
                    error([mfilename ':badnorm'], 'Unrecognized norm choice.');
            end
        end
        switch params.gauge
            case 'zerosum'
                maxentparams = changegauge(maxentparams, 'zerosum');
            case 'original'
            otherwise
                error([mfilename ':badgauge'], 'Unrecognized gauge choice.');
        end
        
        reduced = blockapply(maxentparams.couplings, params.norm, maxentparams);
    otherwise
        error([mfilename ':badmethod'], 'Unrecognized method.');
end

reduced = zeroblockdiag(reduced, 1);

end