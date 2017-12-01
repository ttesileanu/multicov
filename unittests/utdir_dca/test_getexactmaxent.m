function test_getexactmaxent
% TEST_GETEXACTMAXENT Test getexactmaxent.m.

timer = tic;
stats1.type = 'stats';
stats1.alphabets = {'dna', 'multi3'};
stats1.alphawidths = [3 ; 2];
stats1.refseq = struct('seqdb', {'', ''}, 'seqid', {'', ''}, 'map', {1:3, 1:2});
binmap = getbinmap(stats1);
normalize = @(v) v / (sum(v) + 0.01 + 0.98*rand(1, 1));
stats1.freq1 = cell2mat(cellfun(@(idxs) flatten(normalize(0.01 + 0.98*rand(size(idxs)))), ...
    binmap(:), 'uniform', false));
stats1.freq2 = zeroblockdiag(stats1.freq1(:)*stats1.freq1(:)', stats1) + diag(stats1.freq1);
stats1.cmat = stats1.freq2 - stats1.freq1(:)*stats1.freq1(:)';

params = getexactmaxent(stats1, 'verbose', false, 'tol', 1e-6);
stats1res = getmaxentstatsslow(paramsincludegaps(params), 'verbose', false);

utexpect(statscompare(stats1, stats1res, 'tol', 15*eps), 'getexactmaxent diagonal couplings', timer);

timer = tic;
% creating a covariance matrix that is compatible with freq1 and with the
% structure of the problem, and is otherwise random but is "large" (i.e.,
% not generated from the sampling noise from a finite-sized alignment) is
% actually quite hard...
basepattern.alphabets = {'binary', 'rna'};
basepattern.alphawidths = [2; 1];
basepattern.refseq = struct('seqdb', {'', ''}, 'seqid', {'', ''}, 'map', {1:4, 1:2});
binmap = getbinmap(basepattern);
freq1 = cell2mat(cellfun(@(idxs) flatten(normalize(0.01 + 0.98*rand(size(idxs)))), ...
    binmap(:), 'uniform', false));
baserndalign = alngenrandom(1024, freq1, basepattern);

% from this completely random alignment, generate a correlated alignment
rndalign = baserndalign;
rndalign.alphabets = {'binary', 'rna'};
rndalign.alphawidths = [4; 2];
rndalign.refseq = struct('seqdb', {'', ''}, 'seqid', {'', ''}, 'map', {1:4, 1:2});
rndalign.data = [baserndalign.data(:, 1) baserndalign.data(:, 2) baserndalign.data(:, 2) baserndalign.data(:, 1) ...
    baserndalign.data(:, 3) baserndalign.data(:, 3)];

% complete correlations are boring, and a degenerate case which we can't
% handle; make these fuzzy
for i = 1:size(rndalign.data, 2)
    mask = (rand(size(rndalign.data, 1), 1) < 0.4);
    crtcolmasked = rndalign.data(mask, i);
    crtcolmasked = crtcolmasked(randperm(length(crtcolmasked)));
    rndalign.data(mask, i) = crtcolmasked;
end

% finally, get the stats from this alignment
stats2 = getstats(rndalign, 'freq2', true);

% the point of all this was only to get a target stats structure
params2 = getexactmaxent(stats2, 'verbose', false);
stats2res = getmaxentstatsslow(paramsincludegaps(params2), 'verbose', false);

utexpect(statscompare(stats2, stats2res, 'tol', 1e-3), 'getexactmaxent random couplings', timer);

end