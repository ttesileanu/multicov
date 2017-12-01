function test_statspseudocount
% TEST_STATSPSEUDOCOUNT Test statspseudocount.m.

timer = tic;
N = 256;
n = 32;
alignment1 = alngenrandom(N, n, 'protein');
stats0 = getstats(alignment1, 'freq2', true);
alpha1 = 0.3;
stats1 = statspseudocount(stats0, alpha1);

freq1exp = (1 - alpha1)*stats0.freq1 + alpha1/21;
freq2exp = (1 - alpha1)*stats0.freq2 + alpha1/(21*21);
for i = 1:n
    idxs = 20*(i-1) + (1:20);
    freq2exp(idxs, idxs) = diag(freq1exp(idxs));
end

utexpect(all(abs(flatten(stats1.cmat - stats1.freq2 + stats1.freq1(:)*stats1.freq1(:)')) < eps) && ...
    all(abs(stats1.freq1(:) - freq1exp(:)) < eps) && ...
    all(abs(flatten(stats1.freq2 - freq2exp)) < eps), ...
    'statspseudocount single alphabet, uniform background', timer);

% test different alpha for freq1 vs. freq2 & cmat
timer = tic;
alpha1_diff = 0.1;
freq1exp_diff = (1 - alpha1_diff)*stats0.freq1 + alpha1_diff/21;
stats1_diff = statspseudocount(stats0, alpha1, 'fieldalpha', alpha1_diff);
utexpect(all(abs(stats1_diff.freq1(:) - freq1exp_diff(:)) < eps), ...
    'statspseudocount with fieldalpha', timer);

pattern2.alphabets = {'dna' ; 'rna' ; 'protein'};
pattern2.alphawidths = [16 ; 17 ; 5];
alignment2 = alngenrandom(N, pattern2);
stats0 = getstats(alignment2);
stats0freq2 = stats0.cmat + stats0.freq1(:)*stats0.freq1(:)';
alpha2 = 0.6;
stats2 = statspseudocount(stats0, alpha2, 'background', 'database');

expbkgfreq = cell2mat(arrayfun(@(i) repmat(flatten(getdbbkg(pattern2.alphabets{i})), pattern2.alphawidths(i), 1), ...
    (1:length(pattern2.alphabets))', 'uniform', false));
freq1exp = (1 - alpha2)*stats0.freq1(:) + alpha2*expbkgfreq(:);
freq2exp = (1 - alpha2)*stats0freq2 + alpha2*expbkgfreq(:)*expbkgfreq(:)';
binmap = getbinmap(alignment2);
for i = 1:length(binmap)
    idxs = binmap{i};
    freq2exp(idxs, idxs) = diag(freq1exp(idxs));
end
cmatexp = freq2exp - freq1exp(:)*freq1exp(:)';

utexpect(all(abs(stats2.freq1(:) - freq1exp(:)) < eps) && ...
    all(abs(flatten(stats2.cmat - cmatexp)) < eps), ...
    'statspseudocount multiple alphabets, database background for all', timer);

timer = tic;
gapalignment2 = alnincludegaps(alignment2);
gapstats0 = getstats(gapalignment2);
gapstats2 = statspseudocount(gapstats0, alpha2, 'background', 'database');
stats2gap = statsincludegaps(stats2);
utexpect(statscompare(gapstats2, stats2gap, 'tol', 1e-12), 'statspseudocount with gapped alphabets', timer);

alpha3 = 0.1;
stats3 = statspseudocount(stats0, alpha3, 'background', {'database', 'uniform', 'uniform'});

expbkgfreq = [repmat(flatten(getdbbkg('dna')), pattern2.alphawidths(1), 1) ; ...
    ones(4*pattern2.alphawidths(2), 1)/5 ; ones(20*pattern2.alphawidths(3), 1) / 21];
freq1exp = (1 - alpha3)*stats0.freq1(:) + alpha3*expbkgfreq(:);
freq2exp = (1 - alpha3)*stats0freq2 + alpha3*expbkgfreq(:)*expbkgfreq(:)';
binmap = getbinmap(alignment2);
for i = 1:length(binmap)
    idxs = binmap{i};
    freq2exp(idxs, idxs) = diag(freq1exp(idxs));
end
cmatexp = freq2exp - freq1exp(:)*freq1exp(:)';

utexpect(all(abs(stats3.freq1(:) - freq1exp(:)) < eps) && ...
    all(abs(flatten(stats3.cmat - cmatexp)) < eps), ...
    'statspseudocount multiple alphabets, different background choices', timer);

end