function test_getmaxentparams
% TEST_GETMAXENTPARAMS Test getmaxentparams.m.

timer = tic;
N = 128;
n = 10;
alignment1 = alngenrandom(N, n, 'protein');
dca1 = getdca(alignment1);
params1 = getmaxentparams(dca1, 'extended', false);
invC = -inv(dca1.cmat);
Jexp = zeroblockdiag(invC, alignment1);
utexpect(all(isfield(params1, {'type', 'alphabets', 'alphawidths', 'couplings', 'refseq'})) && ...
    strcmp(params1.type, 'maxent') && isequal(params1.alphabets, alignment1.alphabets) && ...
    isequal(params1.alphawidths, alignment1.alphawidths) && isequal(params1.refseq, alignment1.refseq) && ...
    all(abs(Jexp(:) - flatten(params1.couplings - diag(diag(params1.couplings)))) < eps), ...
    'getmaxentparams single alphabet test, not extended', timer);

% check the diagonal trick
timer = tic;
params1_dtrick = getmaxentparams(dca1, 'extended', false, 'diagtrick', true);
dtrick_diff = diag(params1_dtrick.couplings) - diag(params1.couplings);
utexpect(norm(dtrick_diff) > eps, ...
    'getmaxentparams diagonal trick, trivial test', timer);

utexpect(norm(params1_dtrick.couplings - params1.couplings - diag(dtrick_diff)) < eps, ...
    'getmaxentparams diagonal trick only affects diagonal test', timer);

dtrick_diff_exp = diag(invC) - 2*((invC - Jexp)*dca1.freq1);
utexpect(all(abs(dtrick_diff - dtrick_diff_exp) < 1e-12), ...
    'getmaxentparams diagonal trick, better test');

timer = tic;
pattern2.alphabets = {'dna', 'multi7'};
pattern2.alphawidths = [7 ; 13];
alignment2 = alngenrandom(N, pattern2);
dca2 = getdca(alignment2);
params2 = getmaxentparams(dca2, 'fieldcorrection', true);
params2_nofc = getmaxentparams(dca2, 'fieldcorrection', false);
utexpect(any(abs(diag(params2_nofc.couplings) - diag(params2.couplings)) > eps) && ...
    all(size(params2.couplings) == (7*5+13*7)) && ...
    all(abs(flatten(params2.couplings - params2.couplings')) < 1e-14), ...
    'getmaxentparams multi alphabet trivial test, extended, field correction', timer);

utexpect(all(abs(flatten(params2.couplings - diag(diag(params2.couplings)) - zeroblockdiag(params2.couplings, params2))) < eps) && ...
    all(abs(flatten(params1.couplings - diag(diag(params1.couplings)) - zeroblockdiag(params1.couplings, params1))) < eps), ...
    'getmaxentparams diagonal blocks of coupling matrix are diagonal');

timer = tic;
alignment2_wg = alignment2;
alignment2_wg.alphabets = {'gapdna', 'gapmulti7'};
binmap = getbinmap(alignment2_wg);
testres = true;
for i = 1:length(binmap)
    idx = binmap{i}(1);
    if any(abs(params2.couplings(idx, :)) > eps)
        testres = false;
        break;
    end
end
utexpect(testres, 'getmaxentparams check result is in gaps = 0 gauge', timer);

%params2_noext = getmaxentparams(dca2, 'fieldcorrection', true, 'extended', false);
%params2_nofc_noext = getmaxentparams(dca2, 'fieldcorrection', false, 'extended', false);

timer = tic;
dca2_wg = getdca(alignment2_wg);
freq1_wg = dca2_wg.freq1;
J_nd = params2.couplings - diag(diag(params2.couplings));
corr_exp = -J_nd*freq1_wg(:);
tol = 1e-12;
utexpect(all(abs(flatten(diag(params2.couplings) - diag(params2_nofc.couplings)) - 2*corr_exp(:)) < tol), ...
    'getmaxentparams field correction check', timer);

% test that output is the same if starting from gapped statistics
timer = tic;
gapdca1 = statsincludegaps(dca1);
gapdca2 = statsincludegaps(dca2);
gapparams1 = getmaxentparams(gapdca1, 'extended', false);
gapparams2 = getmaxentparams(gapdca2);
utexpect(all(abs(gapparams1.couplings(:) - params1.couplings(:)) < eps) && ...
    all(abs(gapparams2.couplings(:) - params2.couplings(:)) < eps), ...
    'getmaxentparams independence on gappiness of input', timer);

end