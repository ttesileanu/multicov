function test_statscheck
% TEST_STATSCHECK Test statscheck.m.

alphabet = 'rna';
letters = alphagetletters(alphabet);
N = 128;
n = 32;
data = letters(randi(length(letters), N, n));
utexpect(statscheck(getstats(alnmake(data, alphabet))), 'statscheck getstats from random alignment');

mocks.type = 'stats';
mocks.alphabets = {'protein', 'dna'};
mocks.alphawidths = [2 3];
w = 2*20 + 3*4;
mocks.freq1 = zeros(w, 1);
mocks.cmat = zeros(w, w);
mocks.refseq = struct('seqdb', {'generic', 'generic'}, 'seqid', {'', ''}, ...
    'map', {[1 ; 2], [1 ; 2 ; 3]});
utexpect(statscheck(mocks), 'statscheck mock structure');

utexpect(~statscheck(struct('alphabets', 'a', 'alphawidths', 3)), 'statscheck with non-stats structure');

utexpect(~statscheck(alnmake), 'statscheck with character alignment');

end