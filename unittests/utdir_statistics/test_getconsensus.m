function test_getconsensus
% TEST_GETCONSENSUS Test getconsensus.m.

timer = tic;
N = 128;
n = 64;
alignment = alngenrandom(N, n, 'protein');
binalign = alntobin(alignment);
stats = getstats(binalign);
cons_aln = getconsensus(alignment);
cons_bin = getconsensus(binalign);
cons_sts = getconsensus(stats);
letters = alphagetletters('protein', 'nogap');
[~, cons_num] = max(reshape(stats.freq1, 20, n), [], 1);
cons_exp = letters(cons_num);
utexpect(strcmp(cons_aln, cons_bin) && strcmp(cons_aln, cons_sts) && ...
    strcmp(cons_aln, cons_exp), 'getconsensus single alphabet, any input mode', timer);

timer = tic;
pattern.alphabets = {'rna', 'protein', 'dna'};
pattern.alphawidths = [16 ; 10 ; 32];
alignment2 = alngenrandom(N, pattern);
cons2_exp = char(mode(double(alignment2.data)));
cons2 = getconsensus(alignment2, 'gaps', true);
utexpect(strcmp(cons2, cons2_exp), 'getconsensus multi alphabet, with gaps', timer);

gapalignment = alnincludegaps(alignment);
gapalignment2 = alnincludegaps(alignment2);
gapcons = getconsensus(gapalignment);
gapcons2 = getconsensus(gapalignment2, 'gaps', true);
utexpect(strcmp(gapcons, cons_aln) && strcmp(gapcons2, cons2), ...
    'getconsensus with gapped alphabets');

% and now a case with sequence weights
mockalign = alnmake([...
        'AF--' ; ...
        'AG-V' ; ...
        '---V' ; ...
        'A-F-' ; ...
        'CHFY' ; ...
    ], 'protein');
mockalign.seqw = [1 ; 0.1 ; 0.5 ; 0.2 ; 0.1];
mockcons_exp = 'AFFV';
mockcons = getconsensus(mockalign);
utexpect(strcmp(mockcons, mockcons_exp), 'getconsensus with sequence weights');

% test the numeric output
mockcons_num = getconsensus(mockalign, 'indices', true);
cons2_num = getconsensus(alignment2, 'indices', true);
cons2_ng = getconsensus(alignment2);
binmap2 = getbinmap(alignment2);
cons2_ng_exp = blanks(sum(alignment2.alphawidths));
for i = 1:sum(alignment2.alphawidths)
    letters = alphagetletters(alignment2.alphabets{find(cumsum(alignment2.alphawidths) >= i, 1)}, 'nogap');
    cons2_ng_exp(i) = letters(find(binmap2{i} == cons2_num(i), 1));
end
utexpect(isequal(mockcons_num(:), [1 ; 25 ; 45 ; 78]) && ...
    strcmp(cons2_ng, cons2_ng_exp), 'getconsensus with indices as return values');

end