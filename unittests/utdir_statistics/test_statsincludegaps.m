function test_statsincludegaps
% TEST_STATSINCLUDEGAPS Test statsincludegaps.m.

timer = tic;
N = 128;
pattern.alphabets = {'rna', 'protein', 'binary'};
pattern.alphawidths = [5 ; 7 ; 11];
alignment = alngenrandom(N, pattern);
gapalign = alnincludegaps(alignment);
stats = getstats(alignment, 'freq2', true);
statsgap = getstats(gapalign, 'freq2', true);
[gapstats, changed1] = statsincludegaps(stats);
utexpect(statscompare(statsgap, gapstats), 'statsincludegaps stats structure', timer);

freq1gap = statsincludegaps(stats.freq1, stats);
cmatgap = statsincludegaps(stats.cmat, pattern);
utexpect(all(abs(freq1gap(:) - statsgap.freq1(:)) < eps) && ...
    all(abs(cmatgap(:) - statsgap.cmat(:)) < eps), ...
    'statsincludegaps for vector and matrix');

[~, changed2] = statsincludegaps(gapstats);
utexpect(changed1 && ~changed2, 'statsincludegaps ''changed'' output argument');

end