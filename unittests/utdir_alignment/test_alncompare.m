function test_alncompare
% TEST_ALNCOMPARE Test aln_compare.m

utexpect(alncompare(alnmake, alnmake), 'alncompare empty alignments');

aln0.alphabets = {'protein'};
aln0.alphawidths = 4;
aln0.data = ['ACCW' ; 'FFH-' ; '---D'];
aln0.seqw = [0.5 1 1];
aln0.refseq.seqdb = 'fake';
aln0.refseq.seqid = '1';
aln0.refseq.map = [1 5 6 7];
aln0.annotations = {'rat' ; 'mouse' ; 'T-rex'};
aln0.type = 'character';

utexpect(alncompare(aln0, aln0), 'alncompare identical alignments');

aln0a = aln0;
aln0a.extrafield = 'random';
utexpect(alncompare(aln0, aln0a), 'alncompare alignments differing by irrelevant field');

aln0b = aln0;
aln0b.refseq.seq = 'random';
utexpect(alncompare(aln0, aln0b), 'alncompare alignments differing by irrelevant refseq field');

aln0c = aln0;
aln0c.seqw(1) = aln0c.seqw(1) + eps;
utexpect(alncompare(aln0, aln0c), 'alncompare alignments differing by eps in seqw');

aln1 = aln0;
aln1.seqw(1) = 1;
utexpect(~alncompare(aln0, aln1), 'alncompare different seqw');

aln2 = aln0;
aln2.refseq.map(1) = 2;
utexpect(~alncompare(aln0, aln2), 'alncompare different refseq');

aln3 = aln0;
aln3.annotations{2} = 'owl';
utexpect(~alncompare(aln0, aln3), 'alncompare different annotations');

aln4 = aln0;
aln4.data(2, :) = 'FGH-';
utexpect(~alncompare(aln0, aln4), 'alncompare different data');

aln5 = aln0;
aln5.alphawidhts = 3;
aln5.data = ['ACC' ; 'FFH' ; '---'];
utexpect(~alncompare(aln0, aln5), 'alncompare different width');

aln6 = aln0;
aln6.data(4, :) = 'AADD';
utexpect(~alncompare(aln0, aln6), 'alncompare different height');

aln7 = aln0;
aln7.alphabets{1} = 'dna';
utexpect(~alncompare(aln0, aln7), 'alncompare different alphabet');

aln8 = aln0;
aln8.alphabets = {'protein', 'dna'};
aln8.alphawidths = [2 ; 2];
utexpect(~alncompare(aln0, aln8), 'alncompare different number of alphabets');

aln9 = aln0;
aln9.seqw = aln9.seqw';
aln9.annotations = aln9.annotations';
utexpect(alncompare(aln9, aln0), 'alncompare only transposition differences');

end