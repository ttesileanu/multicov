function results = simmaxent(iniseq, maxentparams, nsteps, varargin)
% SIMMAXENT Simulate maximum entropy model.
%   results = SIMMAXENT(iniseq, maxentparams, nsteps) simulates a maximum
%   entropy model with the given couplings starting at the given sequence.
%   The number of steps can be set to 0, in which case the simulation is
%   run until CTRL+C is pressed.
%
%   The results are placed in a structure with the following fields.
%    'acceptances'
%       Average acceptance rates for the recstep steps up to each stored
%       sequence.
%    'energies'
%       The energies of the stored sequences.
%    'nsteps'
%       Number of steps completed. This can be smaller than the number
%       requested if CTRL+C is used.
%    'sequences'
%       An alignment structure holding the sequences stored every rec-step
%       steps after the burnin. Note that this will always have ungapped
%       alphabets, regardless of the structure of maxentparams.
%
%   Options:
%    'burnin' <n>
%       Number of steps used as burnin (no sequences or energies saved).
%       (default: 0)
%    'dispinterval' <x>
%       Time interval in seconds between progress displays.
%       (default: 10)
%    'recstep' <n>
%       Step to use for recording sequences after the burnin.
%       (default: 1)
%    'verbose' <b>
%       Whether to output progress information or not.
%       (default: true)

% Tiberiu Tesileanu (2014)

if nargin == 1 && strcmp(iniseq, 'hascpp')
    results = (exist('cppsimmaxent', 'file') == 3);
    return;
end

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('burnin', 0, @(x) isnumeric(x) && isscalar(x));
parser.addParamValue('dispinterval', 10, @(x) isnumeric(x) && isscalar(x));
parser.addParamValue('recstep', 1, @(x) isnumeric(x) && isscalar(x));
parser.addParamValue('verbose', true, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

% the C++ function needs these as doubles
params.burnin = double(params.burnin);
params.recstep = double(params.recstep);
params.verbose = double(params.verbose);
params.dispinterval = double(params.dispinterval);

% make sure the gaps are included
maxentparams = paramsincludegaps(maxentparams);

% call the C++ function
alphabetletters = cellfun(@(a) alphagetletters(a, 'nogap'), maxentparams.alphabets, 'uniform', false);
cppres = cppsimmaxent(iniseq, alphabetletters, maxentparams.alphawidths, ...
    maxentparams.couplings, double(nsteps), params);

results = cppres;
seqaln = alnmake;
seqaln.alphabets = alphaexcludegaps(maxentparams.alphabets);
seqaln.alphawidths = maxentparams.alphawidths;
seqaln.refseq = maxentparams.refseq;
seqaln.seqw = ones(size(results.sequences, 1), 1);
seqaln.annotations = repmat({}, size(results.sequences, 1), 1);
seqaln.data = results.sequences;
results.sequences = seqaln;

end