function test_findseqmap
% TEST_FINDSEQMAP Test findseqmap.m.

[~, align] = nwalign('GAACC', 'AUGACACGCUU', 'glocal', true, 'alphabet', 'nt');
mapping = findseqmap(align);
utexpect(isequal(mapping, [3 ; 4 ; 6 ; 7 ; 9]), 'findseqmap with nwalign');

mapping2 = findseqmap('GAACC', 'AUGACACGCUU', 'rna');
mapping2a = findseqmap('GAACC', 'AUGACACGCUU', 'nt');
mapping2b = findseqmap('GAACC', 'AUGACACGCUU', 'gaprna');
utexpect(isequal(mapping, mapping2) && isequal(mapping2, mapping2a) && ...
    isequal(mapping2, mapping2b), 'findseqmap with two nucleic sequences');

mapping3 = findseqmap('QRESGC', 'PPFQTESGCR', 'protein');
mapping3a = findseqmap('QRESGC', 'PPFQTESGCR', 'aa');
mapping3b = findseqmap('QRESGC', 'PPFQTESGCR', 'gapprotein');
utexpect(isequal(mapping3, [4 ; 5 ; 6 ; 7 ; 8 ; 9]) && isequal(mapping3, mapping3a) && ...
    isequal(mapping3, mapping3b), 'findseqmap with two protein sequences');

mapping4 = findseqmap('QATWQSRIQERAGACTQHYA', 'QSRIQERAGACS', 'protein');
utexpect(isequal(mapping4, [0 ; 0 ; 0 ; 0 ; (1:12)' ; 0 ; 0 ; 0 ; 0]), ...
    'findseqmap with some non-matches');

end