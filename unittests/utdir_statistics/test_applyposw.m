function test_applyposw
% TEST_APPLYPOSW Test applyposw.m.

timer = tic;
N = 128;
n = 16;
alignment = alngenrandom(N, n, 'protein');
stats = getstats(alignment);
v1 = randn(20*n, 1);
stats_pw = applyposw(stats, @(f, s) v1);
utexpect(all(all(abs(stats_pw.cmat - stats.cmat .* (v1*v1')) < eps)), ...
    'applyposw full amount', timer);

timer = tic;
fraction = 0.3;
stats = getstats(alignment, 'freq2', true);
stats_partial = applyposw(stats, @(f, s) v1, fraction);
freq1exp = (1-fraction)*stats.freq1(:) + fraction*stats_pw.freq1(:);
cmatexp = (1-fraction)*stats.cmat + fraction*stats_pw.cmat;
utexpect(all(all(abs(stats_partial.cmat - stats_partial.freq2 + stats_partial.freq1(:)*stats_partial.freq1(:)') < eps)) && ...
    all(abs(stats_partial.freq1(:) - freq1exp(:)) < eps) && ...
    all(all(abs(stats_partial.cmat - cmatexp) < eps)), ...
    'applyposw partial', timer);

end