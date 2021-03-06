function test_alntobin
% TEST_ALNTOBIN Test alntobin.m.

alignment1 = alnmake(['ACA' ; 'GUA' ; '-A-'], 'rna');
binalign1 = alntobin(alignment1);
bindata1_exp = [...
    1 0 0 0   0 1 0 0   1 0 0 0 ; ...
    0 0 1 0   0 0 0 1   1 0 0 0 ; ...
    0 0 0 0   1 0 0 0   0 0 0 0 ...
];
utexpect(bincheck(binalign1) && isequal(binalign1.data, bindata1_exp), ...
    'alntobin rna');

alignment2 = alnmake(['DF' ; 'YA' ; '-C'], 'protein');
alignment = alnadd(alignment1, alignment2);
binalign2 = alntobin(alignment2);
binalign = alntobin(alignment);
bindata2_exp = [...
    % A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
      0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ...
      0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 ...
      1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ; ...
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ...
      0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ...
];
utexpect(bincheck(binalign2) && bincheck(binalign) && ...
    isequal(binalign.data, [binalign1.data binalign2.data]) && ...
    isequal(binalign2.data, bindata2_exp), 'alntobin multi-alphabet');

end