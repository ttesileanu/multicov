function test_alnscramble
% TEST_ALNSCRAMBLE Test alnscramble.m.

N = 256;
n = 32;
alignment0 = alngenrandom(N, n, 'protein');
alignment = alnadd(alignment0, alignment0);

timer = tic;
freq1 = getfreq1(alignment);
scrambled = alnscramble(alignment);
freq1_scrambled = getfreq1(scrambled);
utexpect(all(abs(freq1 - freq1_scrambled) < eps) && ...
    any(scrambled.data(:) ~= alignment.data(:)), ...
    'alnscramble', timer);

timer = tic;
rng('default');
scrambled2 = alnscramble(alignment);
stats = getstats(scrambled2);
cmatred = blockapply(stats.cmat, @(m) norm(m, 'fro'), stats);
utexpect(all(abs(flatten(cmatred - diag(diag(cmatred)))) < 20/N), ...
    'anlscramble no correlations', timer);

end