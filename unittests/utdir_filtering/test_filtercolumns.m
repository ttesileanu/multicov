function test_filtercolumns
% TEST_FILTERCOLUMNS Test filtercolumns.m.

thresh = 1/21;

alignment = alngenrandom(256, 64, 'protein');
align_clean = filtercolumns(alignment, thresh);

gapfrac0 = mean(alignment.data == '-', 1);
gapfrac = mean(align_clean.data == '-', 1);
utexpect(max(gapfrac) <= thresh && sum(gapfrac0 <= thresh) == length(gapfrac), ...
    'filtercolumns single alphabet');

align_clean_2 = filtercolumns(alignment, thresh, 'edges', true);
gapfrac0 = mean(alignment.data == '-', 1);
sidx = find(gapfrac0 <= thresh, 1, 'first');
eidx = find(gapfrac0 <= thresh, 1, 'last');
if isempty(sidx)
    idxs = [];
else
    idxs = sidx:eidx;
end

utexpect(isequal(align_clean_2.data, alignment.data(:, idxs)), ...
    'filtercolumns single alphabet, ''edges'' true');

rng(7e5);
pattern.alphabets = {'rna', 'dna', 'protein'};
pattern.alphawidths = [3 ; 2 ; 5];
alignment2 = alngenrandom(256, [ones(5*4, 1)*0.94/4 ; ones(5*20, 1)/21], pattern);
align2_clean = filtercolumns(alignment2, thresh);
gapfrac = mean(align2_clean.data == '-', 1);
utexpect(max(gapfrac) <= thresh, 'filtercolumns multi-alphabet');

gapalignment2 = alnincludegaps(alignment2);
gapalign2_clean = filtercolumns(gapalignment2, thresh);
gapalign2_clean_exp = alnincludegaps(align2_clean);
utexpect(isequal(gapalign2_clean, gapalign2_clean_exp), ...
    'filtercolumns with gapped alphabets');

aln3_1 = alnmake(repmat('A', 256, 5), 'protein');
alignment3 = alnadd(aln3_1, alnmake(repmat('-', 256, 5), 'dna'));
align3_clean = filtercolumns(alignment3, thresh);
utexpect(alncompare(align3_clean, aln3_1), 'filtercolumns with disappearing alphabet');

end