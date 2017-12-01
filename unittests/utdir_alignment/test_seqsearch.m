function test_seqsearch
% TEST_SEQSEARCH Test seqsearch.m.

% make a mock alignment to search in
N = 256;
n = 64;
alignment0 = alngenrandom(N, n, 'protein');
alignment0.annotations = arrayfun(@(i) int2str(3*i), 1:N, 'uniform', false);
alignment0.seqw = rand(N, 1);

alignment1 = seqsearch(alignment0, 'idx', 13, 'verbose', false);
align1_exp = alignment0;
align1_exp.data([1 13], :) = align1_exp.data([13 1], :);
align1_exp.seqw([1 13]) = align1_exp.seqw([13 1]);
align1_exp.annotations([1 13]) = align1_exp.annotations([13 1]);
utexpect(isequal(alignment1, align1_exp), 'seqsearch by index');

[~, idx] = seqsearch(alignment0, 'idx', [13 14 15], 'annot', '42', 'swaptofirst', false, 'verbose', false);
utexpect(idx == 14, 'seqsearch by indices + annotation');

[~, idx] = seqsearch(alignment0, 'annot', ['^' alignment0.annotations{17} '$'], 'swaptofirst', false, 'verbose', false);
utexpect(isscalar(idx) && idx == 17, 'seqsearch by annotation, no swap');

% for sequence search we need a little more control
alignment2 = alnmake(['ACGU' ; '--GG' ; 'A-GG' ; 'GGC-'], 'dna');

[~, ~, idx] = seqsearch(alignment2, 'seq', 'GGC', 'verbose', false);
utexpect(isscalar(idx) && idx == 4, 'seqsearch by sequence');

alignment3.alphabets = {'protein', 'rna'};
alignment3.alphawidths = [6 5];
alignment3.data = [...
    '-AWGGHGG-CA' ; ...
    'D-GG-AC-ACC' ; ...
    'WWGYPDGAACC' ; ...
    'W--IIK-UUAU' ; ...
    '--FDGHUCGAU' ; ...
];
alignment3.refseq(1).seqdb = 'generic';
alignment3.refseq(2).seqdb = 'generic';
alignment3.refseq(1).seqid = '';
alignment3.refseq(2).seqid = '';
alignment3.refseq(1).map = 1:6;
alignment3.refseq(2).map = 1:5;
alignment3.seqw = ones(5, 1);
alignment3.annotations = {'' ; '' ; '' ; '' ; ''};
alignment3.type = 'character';
[newalign, idx, oldidx] = seqsearch(alignment3, 'seq', {'WWGYP', 'GAACC'}, 'verbose', false);

gapalignment3 = alnincludegaps(alignment3);
[gapnewalign, gapidx, gapoldidx] = seqsearch(gapalignment3, 'seq', {'WWGYP', 'GAACC'}, 'verbose', false);
newalign_exp = alnincludegaps(newalign);
utexpect(idx == 1 && isscalar(oldidx) && oldidx == 3 && ...
    all(newalign.data(idx, :) == alignment3.data(oldidx, :)) && ...
    all(newalign.data(oldidx, :) == alignment3.data(idx, :)) && ...
    isequal(newalign_exp, gapnewalign) && isequal(idx, gapidx) && isequal(oldidx, gapoldidx), ...
    'seqsearch by sequence, multiple alphabets, gapped or ungapped alphabets');

end
