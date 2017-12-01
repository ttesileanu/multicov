function test_alnincludegaps
% TEST_ALNINCLUDEGAPS Test alnincludegaps.m.

pattern.alphabets = {'protein', 'gapdna', 'rna'};
pattern.testfield = 'foo';
[gapalign, changed1] = alnincludegaps(pattern);
utexpect(all(isfield(gapalign, {'alphabets', 'testfield'})) && ...
    strcmp(pattern.testfield, 'foo') && ...
    isequal(gapalign.alphabets(:), {'gapprotein' ; 'gapdna' ; 'gaprna'}), ...
    'alnincludegaps alignment-like structure');

timer = tic;
N = 256;
pattern.alphawidths = [5 ; 11; 7];
alignment = alngenrandom(N, pattern);
gapalign = alnincludegaps(alignment);

gapalign_exp = alignment;
gapalign_exp.alphabets = {'gapprotein', 'gapdna', 'gaprna'};
utexpect(isequal(gapalign_exp, gapalign), 'alncinludegaps character alignment', timer);

[~, changed2] = alnincludegaps(gapalign);
utexpect(changed1 && ~changed2, 'alnincludegaps ''changed'' output argument');

timer = tic;
binalign = alntobin(alignment);
bingapalign = alntobin(gapalign);
gapbinalign = alnincludegaps(binalign);
utexpect(isequal(bingapalign.alphabets, gapbinalign.alphabets) && ...
    isequal(bingapalign.data, gapbinalign.data), ...
    'alnincludegaps binary alignment', timer);

end