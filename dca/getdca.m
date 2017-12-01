function stats = getdca(alignment, varargin)
% GETDCA Calculate DCA matrix.
%   stats = GETDCA(alignment) calculates the DCA matrix, and returns it
%   along with single-site statistics for the alignment. This is just a
%   usual statistics structures, as returned by getstats, supplemented by a
%   statistical regularization procedure.
%
%   stats = GETDCA(stats0) uses a pre-calculated statistics structure.
%
%   Options:
%    'freq2' <b>
%       Whether to keep the 'freq2' information or not.
%       (default: no, to save space)
%    'regalpha' <x>
%       Value to use for the regularization parameter. Set to zero to not
%       perform the regularization.
%       (default: 0.5)
%    'regfct' <f>
%       Regularization function to use (a function handle). This is called
%       with a statistics structure and an amount (see 'regalpha') as
%       parameters. Set this to an empty matrix to skip regularization.
%       (default: @statspseudocount)
%    'small_fields_pc' <b>
%       Set to true to use a different (small) pseudocount for the fields
%       compared to the couplings. Note that this passes the arguments
%       'fieldalpha', 1/n_seqs to the 'regfct' function.

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('freq2', false, @(b) isscalar(b) && islogical(b));
parser.addParamValue('regalpha', 0.5, @(x) isscalar(x) && isnumeric(x));
parser.addParamValue('regfct', @statspseudocount, @(x) isempty(x) || isa(x, 'function_handle'));
parser.addParamValue('small_fields_pc', false, @(b) isscalar(b) && islogical(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

if isempty(params.regfct)
    params.regalpha = 0;
end

stats = getstats(alignment, 'freq2', params.freq2);

if params.regalpha > 0
    if params.small_fields_pc
        extra_args = {'fieldalpha', 1/stats.nseqs};
    else
        extra_args = {};
    end
    stats = params.regfct(stats, params.regalpha, extra_args{:});
end

end