function resparams = getmaxentparams(stats, varargin)
% GETMAXENTPARAMS Estimate the fields and couplings in a maximum entropy
% model from data.
%   params = GETMAXENTPARAMS(stats) estimates the couplings for a maximum
%   entropy model that recovers the 1- and 2-point correlation functions
%   from the given alignment. The values of the fields h_i(A) are stored on
%   the diagonal of the returned matrix, such that couplings(a, a) = 2*h_a.
%
%   The estimation is done using the mean-field approximation described in
%   Morcos et al. 2011, which essentially amounts to approximating the
%   coupling matrix by the inverse of the covariance matrix. The fields can
%   be estimated using the independent-sites approximation (in which case
%   they are essentially the logarithms of the residue frequencies), but by
%   default a correction involving the couplings is enabled (see
%   'fieldcorrection' below). The calculations are done in the 'gaps = 0'
%   gauge, in which all the fields and couplings involving gaps are set to
%   zero to fix the symmetry inherent in the problem (as described in
%   Morcos et al. 2011).
%
%   The parameters are returned in a structure that contains the following
%   fields:
%    'type':
%       Set to 'maxent'.
%    'alphabets':
%    'alphawidths':
%    'refseq':
%       Fields identifying the structure of the coupling matrix and the
%       mapping to the reference sequence(s). If 'extended' is false, these
%       are directly copied from the input. If it is true, they are
%       'extended' by adding 'gap' to each alphabet. This essentially has
%       the effect of signalling to subsequent functions using the
%       parameters that they do contain the values for the gap entries.
%    'couplings':
%       A matrix of couplings J_{ij}(A, B). This includes the fields on the
%       diagonal, J_{ii}(A, A) = 2*h_i(A). Note the factor of 2 (this is
%       justified by the way the energy is written; see getmaxentenergies).
%
%   Options:
%    'extended' <b>
%       If this is true, the coupling matrix and the field vector are
%       extended to include entries for gaps. This is signalled in the
%       'alphabets' array by adding 'gap' to the alphabet names. If
%       'extended' is false, the parameters are returned without the
%       entries corresponding to gaps. There is no loss of information, as
%       the gap entries are redundant.
%       (default: true)
%    'fieldcorrection' <b>
%       If this is true, the fields are calculated to first order in the
%       magnitude of the couplings, which adds a correction compared to the
%       independent-site approximation. Set this to false to ignore the
%       correction.
%       (default: true)
%    'diagtrick' <b>
%       If this is true, the 'diagonal trick' is used, in which the _ii
%       components of the inverse of the covariance matrix are kept in the
%       couplings (and included for the field correction).
%       (default: false)
%
% See also: GETDCA, GETSTATS, GETMAXENTENERGIES, INCLUDEGAPS, EXCLUDEGAPS.

% Tiberiu Tesileanu (2012-2015)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('extended', true, @(b) islogical(b) && isscalar(b));
parser.addParamValue('fieldcorrection', true, @(b) islogical(b) && isscalar(b));
parser.addParamValue('fcfactor', 1, @(x) isscalar(x) && isnumeric(x));
parser.addParamValue('diagtrick', false, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

if ~statscheck(stats)
    error([mfilename ':badarg'], 'First argument should be a statistics structure.');
end

% make sure we don't start with gappy stats
stats = statsexcludegaps(stats);

Pi = stats.freq1;
Cij = stats.cmat;

% the couplings are simply equal to minus the inverse of the covariance
% matrix, deleting the diagonal elements unless we're using the diagonal trick
couplings = -inv(Cij);
diaginv = diag(couplings);
% if we don't use the diagonal trick, we need to get rid of the
% block-diagonal elements of the couplings first
if ~params.diagtrick
    couplings = zeroblockdiag(couplings, stats);
end

% for the fields, I need to normalize the probabilities by P_i(q) for each
% residue...
alnmap = getbinmap(stats);
fields = zeros(size(Pi));
for i = 1:length(alnmap)
    idxs = alnmap{i};
    subPi = Pi(idxs);
    % max to avoid dividing by zero
    PiGap = max(eps, 1 - sum(subPi));
    
    fields(idxs) = log(subPi / PiGap);
end

% do we want the correction?
if params.fieldcorrection
    fields = fields - params.fcfactor*couplings*Pi;
end

if params.diagtrick
    % otherwise we remove the block-diagonal elements now
    couplings = zeroblockdiag(couplings, stats);
    % and add in the diagonal elements of the couplings
    fields = fields + 0.5*diaginv;
end

% put the fields on the diagonal of the coupling matrix
couplings = couplings + 2*diag(fields);

resparams.type = 'maxent';
resparams.couplings = couplings;
resparams.alphabets = stats.alphabets;
resparams.alphawidths = stats.alphawidths;
resparams.refseq = stats.refseq;

if params.extended
    resparams = paramsincludegaps(resparams);
end

end