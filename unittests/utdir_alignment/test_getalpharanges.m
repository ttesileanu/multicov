function test_getalpharanges
% TEST_GETALPHARANGES Test getalpharanges.m.

ranges1 = getalpharanges([5 3 8]);
ranges1_exp = [1 6 9 ; 5 8 16];
utexpect(isequal(ranges1, ranges1_exp), 'getalpharanges with alphawidths');

alignment2 = alnmake(['ACGU' ; '--GG' ; 'A-GG' ; 'GGC-'], 'dna');
utexpect(isequal(getalpharanges(alignment2), [1 ; 4]), 'getalpharanges with alignment');

% test for binary ranges
alignment1 = alnmake(['AAWYML' ; 'DEGHPQ' ; 'CTLHKW' ; 'A--FG-'], 'protein');
alignment = alnadd(alignment1, alignment2);
utexpect(isequal(getalpharanges(alignment), [1 7 ; 6 10]) && ...
    isequal(getalpharanges(alignment, 'binary'), [1 121 ; 120 136]), ...
    'getalpharanges returning binary ranges');

end