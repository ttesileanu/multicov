function test_simmaxent
% TEST_SIMMAXENT Test simmaxent.m.

global multicov_ut_cpppolicy;

if ~simmaxent('hascpp')
    if strcmp(multicov_ut_cpppolicy, 'warn')
        warning([mfilename ':nocpp'], 'Cannot find C++ extension for simmaxent. Simmaxent requires this extension to run; please compile the extension.');
    elseif strcmp(multicov_ut_cpppolicy, 'fail')
        utexpect(false, 'simmaxent C++ extension missing, simmaxent cannot run');
    end
    return;
end

tol = 1e-12;

% invent some couplings
timer = tic;
params.type = 'maxent';
params.alphabets = {'gapprotein', 'gapdna', 'gapmulti5'};
params.alphawidths = [5 ; 11 ; 8];
params.refseq = struct; % don't really care about the contents of refseq here

binmap = getbinmap(params);
len = binmap{end}(end);
m = randn(len, len);
m = m + m';
params.couplings = zeroblockdiag(m, params) + diag(diag(m));

% get a random starting sequence
rng('default');
inialn = alngenrandom(1, params);

% run the Monte Carlo!
nsteps = 1000;
rng('default');
mcres1 = simmaxent(inialn.data, params, nsteps);

% the sizes of all the outputs should agree
utexpect(length(mcres1.acceptances) == nsteps && ...
    length(mcres1.energies) == nsteps && ...
    size(mcres1.sequences.data, 1) == nsteps, ...
    'simmaxent output structure fields have consistent sizes');

% the number of steps should be right
utexpect(abs(nsteps - mcres1.nsteps) < eps, 'simmaxent correct nsteps');

% the acceptances should match the times when the sequence changed
acceptances_exp = false(size(mcres1.acceptances));
for i = 2:length(acceptances_exp)
    acceptances_exp(i) = any(mcres1.sequences.data(i, :) ~= mcres1.sequences.data(i-1, :));
end
utexpect(all((acceptances_exp & (abs(1 - mcres1.acceptances) < eps)) | ...
    (~acceptances_exp & (abs(mcres1.acceptances) < eps))), ...
    'simmaxent matches sequences when recstep = 1');

% test alphabets are ungapped
[~, changed] = alphaexcludegaps(mcres1.sequences);
utexpect(~changed, 'simmaxent sequences field has ungapped alphabets');

% the sequences should be in the correct alphabets
alpharanges = getalpharanges(mcres1.sequences);
testres = true;
for a = 1:length(mcres1.sequences.alphabets)
    block = mcres1.sequences.data(:, alpharanges(1, a):alpharanges(2, a));
    letters = alphagetletters(mcres1.sequences.alphabets{a}, 'nogap');
    if any(~ismember(block, letters))
        testres = false;
        break;
    end
end
utexpect(testres, 'simmaxent sequences are in correct alphabets');

% energies should match what getmaxentenergies gives
energies_exp = getmaxentenergies(mcres1.sequences, params);
utexpect(all(abs(mcres1.energies(:) - energies_exp(:)) < tol), ...
    'simmaxent sequence energies match getmaxentenergies', timer);

% check gauge invariances
timer = tic;
params_a = changegauge(params, 'zerosum');
params_b = changegauge(params, 'gap');
rng('default');
mcres1_a = simmaxent(inialn.data, params_a, nsteps);
rng('default');
mcres1_b = simmaxent(inialn.data, params_b, nsteps);
utexpect(all(mcres1.sequences.data(:) == mcres1_a.sequences.data(:)) && ...
    all(mcres1.sequences.data(:) == mcres1_b.sequences.data(:)), ...
    'simmaxent gauge invariance', timer);

% test burnin and recstep
timer = tic;
rng('default');
mcres2 = simmaxent(inialn.data, params, nsteps, 'burnin', 500, 'recstep', 5);

utexpect(length(mcres2.acceptances) == length(mcres2.energies) && ...
    length(mcres2.energies) == size(mcres2.sequences.data, 1) && ...
    length(mcres2.energies) == 100, ...
    'simmaxent result with burnin&recstep has consistent field sizes');

utexpect(all(all(mcres2.sequences.data == mcres1.sequences.data(501:5:end, :))), ...
    'simmaxent with burnin&recstep sequence storage OK', timer);

timer = tic;
energies2_exp = getmaxentenergies(mcres2.sequences, params);
utexpect(all(abs(mcres2.energies(:) - energies2_exp(:)) < tol), ...
    'simmaxent with burnin&recstep sequence energies match getmaxentenergies', timer);

% the acceptances should be between 0 and 1
utexpect(all(mcres2.acceptances >= 0 & mcres2.acceptances <= 1), ...
    'simmaxent acceptances in correct range for recstep>1');

% make a simple test with uncoupled sites
timer = tic;

normalize_freq = @(v) v/sum(v);
freq1 = cell2mat(cellfun(@(idxs) normalize_freq(0.02 + 0.96*rand(length(idxs), 1)), ...
    binmap(:), 'uniform', false));
params3 = params;
params3.couplings = diag(2*log(freq1));

rng('default');
mcres3 = simmaxent(inialn.data, params3, 1e6, 'burnin', 1e5, 'recstep', 50);
stats3 = getstats(alnincludegaps(mcres3.sequences));

% the tolerance we expect here is at most of order 1 / sqrt(# of sequences)
mctol = 1/sqrt(size(mcres3.sequences, 1));

utexpect(all(abs(stats3.freq1(:) - freq1(:)) < mctol) && ...
    all(flatten(abs(zeroblockdiag(stats3.cmat, stats3))) < mctol), ...
    'simmaxent simulation with uncoupled sites', timer);

end