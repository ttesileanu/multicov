function test_includegaps
% TEST_INCLUDEGAPS Test includegaps.m.

[res1, changed1] = includegaps('protein');
utexpect(strcmp(res1, 'gapprotein'), 'includegaps with single alphabet');

pattern.alphabets = {'dna', 'gapbinary'};
pattern.alphawidths = [16 ; 7];

[res2, changed2] = includegaps(pattern.alphabets);
utexpect(all(strcmp(res2(:), flatten(alphaincludegaps(pattern.alphabets)))), ...
    'includegaps with cell array of alphabets');

[res3, changed3] = includegaps(pattern);
utexpect(isequal(res3, alphaincludegaps(pattern)), ...
    'includegaps with type-less structure');

timer = tic;
N = 64;
alignment = alngenrandom(N, pattern);
[res4, changed4] = includegaps(alignment);
utexpect(isequal(res4, alnincludegaps(alignment)), ...
    'includegaps with character alignment', timer);

timer = tic;
binalign = alntobin(alignment);
[res5, changed5] = includegaps(binalign);
utexpect(isequal(res5, alnincludegaps(binalign)), ...
    'includegaps with binary alignment', timer);

timer = tic;
stats = getstats(alignment);
[res6, changed6] = includegaps(stats);
utexpect(isequal(res6, statsincludegaps(stats)), ...
    'includegaps with statistics structure', timer);

timer = tic;
alignment2 = alngenrandom(N, 9, 'protein');
stats2 = getstats(alignment2);
dca2 = getdca(stats2);
params2 = getmaxentparams(dca2, 'extended', false);
[res7, changed7] = includegaps(params2);
utexpect(isequal(res7, paramsincludegaps(params2)), ...
    'includegaps with maxent params', timer);

[~, changed1a] = includegaps(res1);
[~, changed2a] = includegaps(res2);
[~, changed3a] = includegaps(res3);
[~, changed4a] = includegaps(res4);
[~, changed5a] = includegaps(res5);
[~, changed6a] = includegaps(res6);
[~, changed7a] = includegaps(res7);

utexpect(...
    all([changed1 changed2 changed3 changed4 changed5 changed6 changed7]) && ...
    ~any([changed1a changed2a changed3a changed4a changed5a changed6a changed7a]), ...
    'includegaps ''changed'' output argument');

end