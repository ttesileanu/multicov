function test_getdca_generic
% TEST_GETDCA_GENERIC Some generic tests for getdca.m.

timer = tic;
N = 384;
n = 32;
alignment1 = alngenrandom(N, n, 'protein');
binalign1 = alntobin(alignment1);
stats1 = getstats(binalign1);
dca1 = getdca(alignment1);
dca1a = getdca(binalign1);
dca1b = getdca(stats1);
utexpect(statscompare(dca1, dca1a) && statscompare(dca1, dca1b) && ~statscompare(stats1, dca1), ...
    'getdca same from character, binary alignment, or stats structure', timer);

timer = tic;
dca1_wf2 = getdca(alignment1, 'freq2', true);
utexpect(isfield(dca1_wf2, 'freq2') && ...
    all(abs(flatten(dca1_wf2.cmat + dca1_wf2.freq1(:)*dca1_wf2.freq1(:)' - dca1_wf2.freq2)) < eps), ...
    'getdca ''freq2'' true', timer);

dca1exp = statspseudocount(stats1, 0.5);
utexpect(statscompare(dca1, dca1exp), 'getdca default options');

% small pseudocount for fields
dca1_sfpc = getdca(alignment1, 'small_fields_pc', true);
dca1exp_sfpc = statspseudocount(stats1, 0.5, 'fieldalpha', 1/stats1.nseqs);
utexpect(statscompare(dca1_sfpc, dca1exp_sfpc), ...
    'getdca with small pseudocount for fields');

timer = tic;
pattern2.alphabets = {'rna', 'protein', 'multi5'};
pattern2.alphawidths = [13 ; 5 ; 11];
alignment2 = alngenrandom(N, pattern2);
stats2 = getstats(alignment2);
dca2 = getdca(stats2, 'regalpha', 0.3, 'regfct', @statsshrink);
dca2exp = statsshrink(stats2, 0.3);
utexpect(statscompare(dca2, dca2exp), 'getdca multi alphabet, different ''regalpha'' and ''regfct''', timer);

end
