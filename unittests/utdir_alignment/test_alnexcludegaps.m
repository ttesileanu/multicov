function test_alnexcludegaps
% TEST_ALNEXCLUDEGAPS Test alnexcludegaps.m.

pattern.alphabets = {'gapprotein', 'dna', 'gaprna'};
pattern.testfield = 'foo';
[ungapalign, changed1] = alnexcludegaps(pattern);
utexpect(all(isfield(ungapalign, {'alphabets', 'testfield'})) && ...
    strcmp(pattern.testfield, 'foo') && ...
    isequal(ungapalign.alphabets(:), {'protein' ; 'dna' ; 'rna'}), ...
    'alnexcludegaps character alignment');

timer = tic;
N = 256;
pattern.alphawidths = [5 ; 11; 7];
alignment = alngenrandom(N, pattern);
ungapalign = alnexcludegaps(alignment);
ungapalign_exp = alignment;
ungapalign_exp.alphabets = {'protein', 'dna', 'rna'};
utexpect(isequal(ungapalign_exp, ungapalign), 'alnexcludegaps character alignment', timer);

[~, changed2] = alnexcludegaps(ungapalign);
utexpect(changed1 && ~changed2, 'alnexcludegaps ''changed'' output argument');

timer = tic;
binalign = alntobin(alignment);
binungapalign = alntobin(ungapalign);
ungapbinalign = alnexcludegaps(binalign);
utexpect(isequal(binungapalign.alphabets, ungapbinalign.alphabets) && ...
    isequal(binungapalign.data, ungapbinalign.data), ...
    'alnexcludegaps binary alignment', timer);

end