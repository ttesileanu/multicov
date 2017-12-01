function test_selectcol
% TEST_SELECTCOL Test selectcol.m.

% test for character alignment
N = 256;
pattern.alphabets = {'rna', 'multi7', 'multi13', 'dna'};
pattern.alphawidths = [7 ; 3 ; 10 ; 5];
alignment = alngenrandom(N, pattern);
idxs = [2:4 6 10 9 8 12:2:20];
subalignment = selectcol(alignment, idxs);

subalignment_exp = alignment;
subalignment_exp.alphabets = alignment.alphabets(1:3);
subalignment_exp.alphawidths = [4 ; 3 ; 5];
subalignment_exp.data = alignment.data(:, idxs);
subalignment_exp.refseq = alignment.refseq(1:3);
subalignment_exp.refseq(1).map = alignment.refseq(1).map([2:4 6]);
subalignment_exp.refseq(2).map = alignment.refseq(2).map([3 2 1]);
subalignment_exp.refseq(3).map = alignment.refseq(3).map(2:2:10);

utexpect(alncompare(subalignment, subalignment_exp), ...
    'selectcol for alignment');

% test for binary alignment
binalign = alntobin(alignment);
subbinalign = selectcol(binalign, idxs);
subbinalign_exp = alntobin(subalignment);

utexpect(isequal(subbinalign.data, subbinalign_exp.data), ...
    'selectcol for binary alignment');

% test for statistics structure
stats = getstats(binalign);
substats = selectcol(stats, idxs);
substats_exp = getstats(subbinalign);

utexpect(statscompare(substats, substats_exp), ...
    'selectcol for statistics structure');

% when we have freq2
stats = getstats(binalign, 'freq2', true);
substats = selectcol(stats, idxs);
substats_exp = getstats(subbinalign, 'freq2', true);

utexpect(statscompare(substats, substats_exp), ...
    'selectcol for statistics structure with ''freq2''');

% test with logical selection
mask = [true true false false true false true false true false ...
    false false false false false false false false false false ...
    true false true true false];
subalignment2 = selectcol(alignment, mask);

subalignment2_exp = alignment;
subalignment2_exp.alphabets = alignment.alphabets([1 2 4]);
subalignment2_exp.alphawidths = [4 ; 1 ; 3];
subalignment2_exp.data = alignment.data(:, mask);
subalignment2_exp.refseq = alignment.refseq([1 2 4]);
subalignment2_exp.refseq(1).map = alignment.refseq(1).map(mask(1:7));
subalignment2_exp.refseq(2).map = alignment.refseq(2).map(mask(8:10));
subalignment2_exp.refseq(3).map = alignment.refseq(3).map(mask(21:25));

utexpect(alncompare(subalignment2, subalignment2_exp), ...
    'selectcol for alignment using logical mask');

end