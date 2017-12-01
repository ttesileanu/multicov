function test_alncheck
% TEST_ALNCHECK Test alncheck.m

utexpect(alncheck(alnmake), 'alncheck empty alnmake');

alphabet = 'rna';
letters = alphagetletters(alphabet);
N = 128;
n = 32;
data = letters(randi(length(letters), N, n));
utexpect(alncheck(alnmake(data, alphabet)), 'alncheck alnmake with data');

utexpect(~alncheck(struct('alphabets', 'a', 'alphawidths', 3)), 'alncheck with non-alignment');

utexpect(~alncheck(alntobin(alnmake)), 'alncheck with binary alignment');

end