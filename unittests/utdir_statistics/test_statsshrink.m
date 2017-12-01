function test_statsshrink
% TEST_STATSSHRINK Test statsshrink.m.

timer = tic;
N = 256;
n = 32;
alignment1 = alngenrandom(N, n, 'protein');
stats0 = getstats(alignment1, 'freq2', true);
alpha1 = 0.2;
stats1 = statsshrink(stats0, alpha1);

cmatbkg = zeros(20*n, 20*n);
subfreq = ones(20, 1)/21;
for i = 1:n
    idxs = 20*(i-1) + (1:20);
    cmatbkg(idxs, idxs) = diag(subfreq) - subfreq*subfreq';
end

freq1exp = (1 - alpha1)*stats0.freq1 + alpha1/21;
cmatexp = alpha1*cmatbkg + (1 - alpha1)*stats0.cmat;

utexpect(all(abs(flatten(stats1.cmat - cmatexp)) < eps) && ...
    all(abs(stats1.freq1(:) - freq1exp(:)) < eps) && ...
    all(abs(flatten(stats1.freq2 - cmatexp - stats1.freq1(:)*stats1.freq1(:)')) < eps), ...
    'statsshrink single alphabet, uniform background', timer);

pattern2.alphabets = {'dna' ; 'rna' ; 'protein'};
pattern2.alphawidths = [16 ; 17 ; 5];
alignment2 = alngenrandom(N, pattern2);
stats0 = getstats(alignment2);
alpha2 = 0.5;
stats2 = statsshrink(stats0, alpha2, 'background', 'database');

expbkgfreq = cell2mat(arrayfun(@(i) repmat(flatten(getdbbkg(pattern2.alphabets{i})), pattern2.alphawidths(i), 1), ...
    (1:length(pattern2.alphabets))', 'uniform', false));
freq1exp = (1 - alpha2)*stats0.freq1(:) + alpha2*expbkgfreq(:);
binmap = getbinmap(alignment2);
cmatbkg = zeros(binmap{end}(end));
for i = 1:length(binmap)
    idxs = binmap{i};
    subfreq = flatten(expbkgfreq(idxs));
    cmatbkg(idxs, idxs) = diag(subfreq) - subfreq*subfreq';
end
cmatexp = alpha2*cmatbkg + (1 - alpha2)*stats0.cmat;

utexpect(all(abs(stats2.freq1(:) - freq1exp(:)) < eps) && ...
    all(abs(flatten(stats2.cmat - cmatexp)) < eps), ...
    'statsshrink multiple alphabets, database background for all', timer);

timer = tic;
gapalignment2 = alnincludegaps(alignment2);
gapstats0 = getstats(gapalignment2);
gapstats2 = statsshrink(gapstats0, alpha2, 'background', 'database');
stats2gap = statsincludegaps(stats2);
utexpect(statscompare(gapstats2, stats2gap, 'tol', 1e-12), 'statsshrink with gapped alphabets', timer);

alpha3 = 0.75;
stats3 = statsshrink(stats0, alpha3, 'background', {'database', 'uniform', 'uniform'});

expbkgfreq = [repmat(flatten(getdbbkg('dna')), pattern2.alphawidths(1), 1) ; ...
    ones(4*pattern2.alphawidths(2), 1)/5 ; ones(20*pattern2.alphawidths(3), 1) / 21];
freq1exp = (1 - alpha3)*stats0.freq1(:) + alpha3*expbkgfreq(:);
binmap = getbinmap(alignment2);
cmatbkg = zeros(binmap{end}(end));
for i = 1:length(binmap)
    idxs = binmap{i};
    subfreq = flatten(expbkgfreq(idxs));
    cmatbkg(idxs, idxs) = diag(subfreq) - subfreq*subfreq';
end
cmatexp = alpha3*cmatbkg + (1 - alpha3)*stats0.cmat;

utexpect(all(abs(stats3.freq1(:) - freq1exp(:)) < eps) && ...
    all(abs(flatten(stats3.cmat - cmatexp)) < eps), ...
    'statsshrink multiple alphabets, different background choices', timer);

end