function test_filterrows
% TEST_FILTERROWS Test filterrows.m.

thresh = 1/21;

alignment = alngenrandom(256, 64, 'protein');
align_clean = filterrows(alignment, 'gaps', thresh);

gapfrac0 = mean(alignment.data == '-', 2);
gapfrac = mean(align_clean.data == '-', 2);
utexpect(max(gapfrac) <= thresh && sum(gapfrac0 <= thresh) == length(gapfrac), ...
    'filterrows gap threshold');

alignment2 = alnmake([...
    'AAAAAAAAAA' ; ...
    'ACAAAAAAAA' ; ...
    'AAA-DAAAAA' ; ...
    'AADAAEAYAA' ; ...
    'A--AFGAAWY'], 'protein');
align2_clean = filterrows(alignment2, 'minrefsim', 0.85);
align2_clean_exp = alnmake([...
    'AAAAAAAAAA' ; ...
    'ACAAAAAAAA'], 'protein');
utexpect(alncompare(align2_clean, align2_clean_exp), 'filterrows refdist threshold');

gapalignment2 = alnincludegaps(alignment2);
gapalign2_clean = filterrows(gapalignment2, 'minrefsim', 0.85);
gapalign2_clean_exp = alnincludegaps(align2_clean);
utexpect(isequal(gapalign2_clean, gapalign2_clean_exp), 'filterrows gapped alphabets');

align2_clean2 = filterrows(alignment2, 'gaps', 0.05, 'minrefsim', 0.75);
align2_clean2_exp = alnmake([...
    'AAAAAAAAAA' ; ...
    'ACAAAAAAAA'], 'protein');
utexpect(alncompare(align2_clean2, align2_clean2_exp), 'filterrows gaps&refdist threshold');

end