function test_getstats
% TEST_GETSTATS Test getstats.m.

% check that it contains the right fields
timer = tic;
N = 128;
n = 32;
alignment1 = alngenrandom(N, n, 'protein');
freq1 = getfreq1(alignment1);
freq2 = getfreq2(alignment1);
stats = getstats(alignment1);
utexpect(statscheck(stats) && all(isfield(stats, {'type', 'freq1', 'cmat', 'alphabets', 'alphawidths'})) && ...
    isequal(stats.alphabets, alignment1.alphabets) && isequal(stats.alphawidths, alignment1.alphawidths) && ...
    ~isfield(stats, 'freq2') && strcmp(stats.type, 'stats') && ...
    all(abs(stats.freq1(:) - freq1(:)) < eps) && ...
    all(all(abs(stats.cmat + freq1(:)*freq1(:)' - freq2) < eps)), ...
    'getstats single alphabet, basic test', timer);

stats_wf2 = getstats(stats, 'freq2', true);
utexpect(isfield(stats_wf2, 'freq2') && all(abs(freq2(:) - stats_wf2.freq2(:)) < eps), ...
    'getstats add freq2 field');

% check with multiple alphabets and returning freq2
timer = tic;
pattern.alphabets = {'dna', 'rna'};
pattern.alphawidths = [5 ; 4];
alignment2 = alngenrandom(N, pattern);
stats = getstats(alignment2, 'freq2', true);
binalign2 = alntobin(alignment2);
covexp = cov(binalign2.data, 1);
utexpect(isfield(stats, 'freq2') && all(all(abs(stats.cmat + stats.freq1(:)*stats.freq1(:)' - stats.freq2) < eps)) && ...
    all(abs(stats.cmat(:) - covexp(:)) < eps), 'getstats multi alphabet test, with freq2', timer);

stats_nf2 = getstats(stats, 'freq2', false);
utexpect(~isfield(stats_nf2, 'freq2'), 'getstats remove freq2 field');

utexpect(statscompare(stats, getstats(stats, 'freq2', true)), ...
    'getstats keep freq2 if ''freq2'' is true');

end