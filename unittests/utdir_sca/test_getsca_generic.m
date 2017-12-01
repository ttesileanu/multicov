function test_getsca_generic
% TEST_GETSCA_GENERIC Some generic tests for getsca.m.

% some tests for the reduction method
timer = tic;
N = 128;
n = 16;
alignment = alngenrandom(N, n, 'protein');
binalign = alntobin(alignment);
stats = getstats(alignment, 'freq2', true);
sca1 = getsca(alignment);
sca2 = getsca(binalign);
sca3 = getsca(stats);
stscomp = @(s1, s2) all(abs(s1.freq1(:) - s2.freq1(:)) < eps) && ...
    (~isfield(s1, 'freq2') || ~isfield(s2, 'freq2') || all(abs(s1.freq2(:) - s2.freq2(:)) < eps)) && ...
    all(abs(s1.cmat(:) - s2.cmat(:)) < eps) && ...
    (~isfield(s1, 'cmatfull') || ~isfield(s2, 'cmatfull') || all(abs(s1.cmatfull(:) - s2.cmatfull(:)) < eps));
utexpect(stscomp(sca1, sca2) && stscomp(sca1, sca3) && stscomp(sca2, sca3), ...
    'getsca with binary or character alignment, or stats structure', timer);

timer = tic;
sca1_noabs = getsca(alignment, 'abs', false, 'freq2', true, 'full', true);
sca1exp = applyposw(stats, @pwfctdkl);
sca1exp.cmatfull = sca1exp.cmat;
sca1exp.cmat = blockapply(sca1exp.cmat, @(m) norm(m, 'fro'), alignment);
utexpect(isfield(sca1_noabs, 'freq2') && isfield(sca1exp, 'freq2') && ...
    isfield(sca1_noabs, 'cmatfull') && ...
    stscomp(sca1_noabs, sca1exp), 'getsca single alphabet, no abs, with freq2 and cmatfull', timer);

timer = tic;
pattern.alphabets = {'rna', 'protein'};
pattern.alphawidths = [16 ; 6];
alignment2 = alngenrandom(N, pattern);
stats2 = getstats(alignment2);
sca = getsca(alignment2, 'pwfct', @(f, s) log(f + eps), 'redfct', []);
posw = log(stats2.freq1 + eps);
sca_exp.freq1 = stats2.freq1 .* posw;
sca_exp.cmat = abs(stats2.cmat .* (posw(:)*posw(:)'));
utexpect(stscomp(sca, sca_exp), 'getsca multi alphabet, pwfct, no redfct', timer);

sca0 = getsca(stats2, 'pwfct', [], 'redfct', [], 'abs', false);
utexpect(stscomp(sca0, stats2), 'getsca with no pwfct, no redfct');

% brief tests for the binary and projection approximations

% get the binary approximation SCA
timer = tic;
sca1bin = getsca(alignment, 'method', 'binary');
sca1bin_fromstats = getsca(stats, 'method', 'binary');
posw = pwfctdkl(stats.freq1, stats);
binalign.data = bsxfun(@times, binalign.data, posw(:)');
idxs = getconsensus(alignment, 'indices', true);
data1bin = binalign.data(:, idxs);
sca1exp = abs(cov(data1bin, 1));
utexpect(all(abs(sca1bin.cmat(:) - sca1exp(:)) < 16*eps) && ...
    stscomp(sca1bin_fromstats, sca1bin), ...
    'getsca binary approximation', timer);

% get the projection approximation SCA
timer = tic;
binalign2 = alntobin(alignment2);
sca2proj = getsca(binalign2, 'method', 'projection', 'abs', false);
sca2proj_fromstats = getsca(stats2, 'method' ,'projection', 'abs', false);
posw = pwfctdkl(stats2.freq1, stats2);
binalign2.data = bsxfun(@times, binalign2.data, posw(:)');
binmap = getbinmap(alignment2);
freq1 = getfreq1(binalign2);
data2proj = zeros(size(alignment2.data));
for i = 1:sum(alignment2.alphawidths)
    idxs = binmap{i};
    vec = freq1(idxs) / norm(freq1(idxs));
    data2proj(:, i) = binalign2.data(:, idxs)*vec(:);
end
sca2exp = cov(data2proj, 1);

utexpect(all(abs(sca2proj.cmat(:) - sca2exp(:)) < 16*eps) && ...
    stscomp(sca2proj, sca2proj_fromstats), ...
    'getsca projection approximation', timer);

end