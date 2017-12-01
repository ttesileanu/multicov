function test_alphaincludegaps
% TEST_ALPHAINCLUDEGAPS Test alphaincludegaps.m.

utexpect(strcmp(alphaincludegaps('xyz'), 'gapxyz') && ...
    strcmp(alphaincludegaps('protein'), 'gapprotein') && ...
    strcmp(alphaincludegaps('gapprotein'), 'gapprotein'), ...
    'alphaincludegaps with single string');

list1 = {'protein'};
list2 = {'dna', 'rna', 'multi5'};
list3 = {'gapdna', 'protein'};
list1exp = {'gapprotein'};
list2exp = {'gapdna', 'gaprna', 'gapmulti5'};
list3exp = {'gapdna', 'gapprotein'};
utexpect(isequal(list1exp(:), flatten(alphaincludegaps(list1))) && ...
    isequal(list2exp(:), flatten(alphaincludegaps(list2))) && ...
    isequal(list3exp(:), flatten(alphaincludegaps(list3))), ...
    'alphaincludegaps with cell array');

structure.alphabets = {'protein', 'gapdna', 'binary'};
structure.alphawidths = [15 ; 7 ; 3];
[gapstructure, changed] = alphaincludegaps(structure);
utexpect(isequal(gapstructure.alphabets(:), {'gapprotein' ; 'gapdna' ; 'gapbinary'}) && ...
    isequal(gapstructure.alphawidths, structure.alphawidths), 'alphaincludegaps with structure');

[~, changed2] = alphaincludegaps(list2);
[~, changed3] = alphaincludegaps(list3);
[~, changed4] = alphaincludegaps('gapprotein');
utexpect(changed && changed2 && changed3 && ~changed4, ...
    'alphaincludegaps ''changed'' output argument');

end