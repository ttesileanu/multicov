function test_statscompare
% TEST_STATSCOMPARE Test statscompare.m.

timer = tic;
N = 256;
n = 32;
alignment1 = alngenrandom(N, n, 'protein');
stats1 = getstats(alignment1, 'freq2', true);
utexpect(statscompare(stats1, stats1), 'statscompare identical structures', timer);

timer = tic;
alignment2 = alngenrandom(N, n, 'protein');
stats2 = getstats(alignment2, 'freq2', true);
utexpect(~statscompare(stats1, stats2), 'statscompare two alignments, same alignment structure', timer);

timer = tic;
stats2a = getstats(alignment2);
utexpect(~statscompare(stats2, stats2a), 'statscompare existence of freq2', timer);

timer = tic;
pattern3.alphabets = {'dna' ; 'rna'};
pattern3.alphawidths = [15 ; 7];
alignment3 = alngenrandom(2*N, pattern3);
stats3 = getstats(alignment3);
[cres3, details3] = statscompare(stats2a, stats3);
utexpect(~cres3 && ...
    all(isfield(details3, {'dalphabets', 'dalphawidths'})) && ...
    all(~isfield(details3, {'dfreq1', 'dcmat', 'dfreq2'})) && ...
    details3.dalphabets && details3.dalphawidths, ...
    'statscompare two alignments, different alignment structure', timer);

timer = tic;
stats1a = stats1;
stats1a.freq1(5) = stats1a.freq1(5) + 0.05;
stats1a.freq2(3, 4) = stats1a.freq2(3, 4) - 0.03;
stats1a.freq2(4, 3) = stats1a.freq2(3, 4);
stats1.cmat(3, 4) = stats1a.cmat(3, 4) + 0.04;
stats1.cmat(4, 3) = stats1a.cmat(3, 4);
[cres1a, details1a] = statscompare(stats1a, stats1, 'tol', 0.06);
utexpect(~statscompare(stats1a, stats1) && cres1a && ...
    all(isfield(details1a, {'dalphabets', 'dalphawidths', 'dfreq1', 'dcmat', 'dfreq2'})) && ...
    ~details1a.dalphabets && ~details1a.dalphawidths && ...
    all(abs([details1a.dfreq1(1) details1a.dfreq2(1) details1a.dcmat(1)] - [0.05 0.03 0.04]) < eps) && ...
    abs(details1a.dfreq1(2) - mean(abs(stats1a.freq1(:) - stats1.freq1(:)))) < eps && ...
    abs(details1a.dfreq2(2) - mean(abs(stats1a.freq2(:) - stats1.freq2(:)))) < eps && ...
    abs(details1a.dcmat(2) - mean(abs(stats1a.cmat(:) - stats1.cmat(:)))) < eps, ...
    'statscompare with different tolerance', timer);

% test with nans
stats1b = stats1;
stats1b.freq1(2) = nan;
utexpect(~statscompare(stats1b, stats1b), 'statscompare with NaNs');

end