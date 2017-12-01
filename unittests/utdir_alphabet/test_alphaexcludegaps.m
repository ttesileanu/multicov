function test_alphaexcludegaps
% TEST_ALPHAEXCLUDEGAPS Test alphaexcludegaps.m.

utexpect(strcmp(alphaexcludegaps('gapxyz'), 'xyz') && ...
    strcmp(alphaexcludegaps('gapprotein'), 'protein') && ...
    strcmp(alphaexcludegaps('protein'), 'protein'), ...
    'alphaexcludegaps with single string');

list1 = {'gapprotein'};
list2 = {'gapdna', 'gaprna', 'gapmulti5'};
list3 = {'dna', 'gapprotein'};
list1exp = {'protein'};
list2exp = {'dna', 'rna', 'multi5'};
list3exp = {'dna', 'protein'};
utexpect(isequal(list1exp(:), flatten(alphaexcludegaps(list1))) && ...
    isequal(list2exp(:), flatten(alphaexcludegaps(list2))) && ...
    isequal(list3exp(:), flatten(alphaexcludegaps(list3))), ...
    'alphaexcludegaps with cell array');

structure.alphabets = {'gapprotein', 'dna', 'gapbinary'};
structure.alphawidths = [15 ; 7 ; 3];
[gapstructure, changed] = alphaexcludegaps(structure);
utexpect(isequal(gapstructure.alphabets(:), {'protein' ; 'dna' ; 'binary'}) && ...
    isequal(gapstructure.alphawidths, structure.alphawidths), 'alphaexcludegaps with structure');

[~, changed2] = alphaexcludegaps(list2);
[~, changed3] = alphaexcludegaps(list3);
[~, changed4] = alphaexcludegaps('binary');
utexpect(changed && changed2 && changed3 && ~changed4, ...
    'alphaexcludegaps ''changed'' output argument');

end