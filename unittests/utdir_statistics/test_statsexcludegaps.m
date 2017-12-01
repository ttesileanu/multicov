function test_statsexcludegaps
% TEST_STATSEXCLUDEGAPS Test statsexcludegaps.m.

timer = tic;
N = 128;
pattern.alphabets = {'rna', 'protein', 'binary'};
pattern.alphawidths = [5 ; 7 ; 11];
alignment = alngenrandom(N, pattern);
gapalign = alnincludegaps(alignment);
stats = getstats(alignment, 'freq2', true);
statsgap = getstats(gapalign, 'freq2', true);
[ungapstats, changed1] = statsexcludegaps(statsgap);
utexpect(statscompare(stats, ungapstats), 'statsexcludegaps stats structure', timer);

freq1ungap = statsexcludegaps(statsgap.freq1, statsgap);
cmatungap = statsexcludegaps(statsgap.cmat, statsgap);
utexpect(all(abs(freq1ungap(:) - stats.freq1(:)) < eps) && ...
    all(abs(cmatungap(:) - stats.cmat(:)) < eps), ...
    'statsexcludegaps for vector and matrix');

[~, changed2] = statsexcludegaps(ungapstats);
utexpect(changed1 && ~changed2, 'statsexcludegaps ''changed'' output argument');

end