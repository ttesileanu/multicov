function test_excludegaps
% TEST_EXCLUDEGAPS Test excludegaps.m.

[res1, changed1] = excludegaps('gapprotein');
utexpect(strcmp(res1, 'protein'), ...
    'excludegaps with single alphabet');

pattern.alphabets = {'gapdna', 'binary'};
pattern.alphawidths = [16 ; 7];

[res2, changed2] = excludegaps(pattern.alphabets);
utexpect(all(strcmp(res2(:), flatten(alphaexcludegaps(pattern.alphabets)))), ...
    'excludegaps with cell array of alphabets');

[res3, changed3] = excludegaps(pattern);
utexpect(isequal(res3, alphaexcludegaps(pattern)), ...
    'excludegaps with type-less structure');

timer = tic;
N = 64;
alignment = alngenrandom(N, pattern);
[res4, changed4] = excludegaps(alignment);
utexpect(isequal(res4, alnexcludegaps(alignment)), ...
    'excludegaps with character alignment', timer);

timer = tic;
binalign = alntobin(alignment);
[res5, changed5] = excludegaps(binalign);
utexpect(isequal(res5, alnexcludegaps(binalign)), ...
    'excludegaps with binary alignment', timer);

timer = tic;
stats = getstats(alignment);
[res6, changed6] = excludegaps(stats);
utexpect(isequal(res6, statsexcludegaps(stats)), ...
    'excludegaps with statistics structure', timer);

timer = tic;
alignment2 = alngenrandom(N, 9, 'protein');
stats2 = getstats(alignment2);
dca2 = getdca(stats2);
params2 = getmaxentparams(dca2);
[res7, changed7] = excludegaps(params2);
utexpect(isequal(res7, paramsexcludegaps(params2)), ...
    'excludegaps with maxent params', timer);

[~, changed1a] = excludegaps(res1);
[~, changed2a] = excludegaps(res2);
[~, changed3a] = excludegaps(res3);
[~, changed4a] = excludegaps(res4);
[~, changed5a] = excludegaps(res5);
[~, changed6a] = excludegaps(res6);
[~, changed7a] = excludegaps(res7);

utexpect(...
    all([changed1 changed2 changed3 changed4 changed5 changed6 changed7]) && ...
    ~any([changed1a changed2a changed3a changed4a changed5a changed6a changed7a]), ...
    'excludegaps ''changed'' output argument');

end