function test_getmaxentstatsslow
% TEST_GETMAXENTSTATSSLOW Test getmaxentstatsslow.m.

rng('default');
tol = 10*eps;

timer = tic;
params.type = 'maxent';
params.alphabets = {'gapmulti5', 'gapbinary'};
params.alphawidths = [2 ; 3];
params.refseq = struct('seqdb', {'generic', 'generic'}, 'seqid', {'', ''}, ...
    'map', {1:2, 1:3});
binmap = getbinmap(params);
prob_normalize = @(v) v / sum(v);
freq_cell = cellfun(@(c) prob_normalize(0.01 + 0.98*rand(length(c), 1)), ...
    binmap, 'uniform', false);
params.couplings = 2*diag(cell2mat(cellfun(@(f) log(f/f(1)), freq_cell(:), 'uniform', false)));

freq = cell2mat(cellfun(@(c) flatten(c(2:end)), freq_cell(:), 'uniform', false));
maxentstats = getmaxentstatsslow(params, 'verbose', false);
utexpect(all(abs(freq(:) - maxentstats.freq1(:)) < tol), ...
    'getmaxentstatsslow diagonal couplings, single-site frequencies check', timer);

utexpect(all(abs(flatten(zeroblockdiag(maxentstats.cmat, maxentstats))) < tol), ...
    'getmaxentstatsslow diagonal couplings, correlation check');

timer = tic;
params2.type = 'maxent';
params2.alphabets = {'gapdna', 'gaprna'};
params2.alphawidths= [2 ; 2];
params2.refseq = struct('seqdb', {'generic', 'generic'}, 'seqid', {'', ''}, ...
    'map', {1:2, 1:2});
binmap2 = getbinmap(params2);
L2 = binmap2{end}(end);
make_symmetric = @(m) (m + m')/2;
params2.couplings = (zeroblockdiag(0.5*make_symmetric(randn(L2, L2)), params2) + 2*diag(randn(L2, 1)));

maxentstats2 = getmaxentstatsslow(params2, 'verbose', false);

global multicov_ut_cpppolicy;

if simmaxent('hascpp')
    % comparison to Monte Carlo, to test what happens when we have non-trivial
    % correlations
    inialn2 = alngenrandom(1, params2);
    mcsim2 = simmaxent(inialn2.data, params2, 1e6, 'burnin', 1e5, 'recstep', 50, ...
        'verbose', false);
    mcstats2 = getstats(mcsim2.sequences, 'freq2', true);
    
    utexpect(statscompare(maxentstats2, mcstats2, 'tol', 0.01), ...
        'getmaxentstatsslow random couplings, comparison to MC', timer);
else
    if strcmp(multicov_ut_cpppolicy, 'warn')
        warning([mfilename ':nocpp'], 'Cannot test getmaxetstatsslow against Monte Carlo results because the C++ extension for simmaxent does not exist.');
    elseif strcmp(multicov_ut_cpppolicy, 'fail')
        utexpect(false, 'C++ extension for simmaxent missing, cannot check against Monte Carlo.');
    end
end

if getmaxentstatsslow('hascpp')
    maxentstats_nocpp = getmaxentstatsslow(params, 'verbose', false, 'nocpp', true);
    utexpect(statscompare(maxentstats, maxentstats_nocpp, 'tol', 10*eps), ...
        'getmaxentstatsslow C++ extension test, diagonal couplings');
    
    maxentstats2_nocpp = getmaxentstatsslow(params2, 'verbose', false, 'nocpp', true);
    utexpect(statscompare(maxentstats2, maxentstats2_nocpp, 'tol', 10*eps), ...
        'getmaxentstatsslow C++ extension test, random couplings');
else
    if strcmp(multicov_ut_cpppolicy, 'warn')
        warning([mfilename ':nocpp'], 'C++ extension for getmaxentstatsslow missing, cannot test it. Compile extension for optimum performance.');
    elseif strcmp(multicov_ut_cpppolicy, 'fail')
        utexpect(false, 'getmaxentstatsslow C++ extension missing.');
    end
end

% XXX should somehow test the partition function output

end