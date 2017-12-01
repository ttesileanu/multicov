function test_bincheck
% TEST_BINCHECK Test bincheck.m

utexpect(bincheck(alntobin(alnmake)), 'bincheck alntobin from empty alnmake');

alphabet = 'rna';
letters = alphagetletters(alphabet);
N = 128;
n = 32;
data = letters(randi(length(letters), N, n));
utexpect(bincheck(alntobin(alnmake(data, alphabet))), 'bincheck alntobin from alnmake with data');

mockba.type = 'binary';
mockba.alphabets = {'protein', 'dna'};
mockba.alphawidths = [2 3];
mockba.data = zeros(128, 2*20 + 3*4);
mockba.seqw = ones(128, 1);
mockba.annotations = repmat({''}, 128, 1);
mockba.refseq = struct('seqdb', {'generic', 'generic'}, 'seqid', {'', ''}, ...
    'map', {[1 ; 2], [1 ; 2 ; 3]});
utexpect(bincheck(mockba), 'bincheck mock alignment');

utexpect(~bincheck(struct('alphabets', 'a', 'alphawidths', 3)), 'bincheck with non-alignment');

utexpect(~bincheck(alnmake), 'bincheck with character alignment');

end