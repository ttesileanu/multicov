function test_alnadd
% TEST_ALNADD Test alnadd.m.

% start with some trivial tests

% make a mock alignment
alphabet = 'protein';
letters = alphagetletters(alphabet);
N = 256;
n = 64;
alignment = alnmake(letters(randi(length(letters), N, n)), alphabet);

aln1 = alnadd(alnmake, alnmake);
aln2 = alnadd(alnmake, alignment);
aln3 = alnadd(alignment, alnmake);
utexpect(alncompare(aln1, alnmake) && alncompare(aln2, alignment) && ...
    alncompare(aln3, alignment), 'alnadd with empty alignments');

% make another mock alignment
alphabet2 = 'dna';
letters2 = alphagetletters(alphabet2);
n2 = 16;
alignment2 = alnmake(letters2(randi(length(letters2), N, n2)), alphabet2);

aln_mix = alnadd(alignment, alignment2);
aln_exp.type = 'character';
aln_exp.alphabets = {alphabet, alphabet2};
aln_exp.alphawidths = [n n2];
aln_exp.seqw = ones(N, 1);
aln_exp.annotations = repmat({''}, N, 1);
aln_exp.refseq = struct('seqdb', {'generic', 'generic'}, ...
    'seqid', {'', ''}, 'map', {(1:n)', (1:n2)'});
aln_exp.data = [alignment.data alignment2.data];
utexpect(alncompare(aln_mix, aln_exp), 'alnadd with different alphabets');

utexpect(alncompare(aln_mix, alnadd(alignment, alignment2.data, alphabet2)), ...
    'alnadd with data matrix');

end