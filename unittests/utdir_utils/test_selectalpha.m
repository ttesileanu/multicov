function test_selectalpha
% TEST_SELECTALPHA Test selectalpha.m.

N = 256;
subalign1 = alngenrandom(N, 16, 'protein');
subalign2 = alngenrandom(N, 7, 'dna');
subalign3 = alngenrandom(N, 13, 'binary');
subalign4 = alngenrandom(N, 3, 'multi36');

alignment1 = alnadd(subalign1, alnadd(subalign2, alnadd(subalign3, subalign4)));

% test for character alignment
alphaselection1 = selectalpha(alignment1, 2);
utexpect(alncompare(alphaselection1, subalign2), 'selectalpha select one alphabet of character alignment');

alphaselection2 = selectalpha(alignment1, [2 4]);
utexpect(alncompare(alphaselection2, alnadd(subalign2, subalign4)), ...
    'selectalpha select multiple alphabets of character alignment');

% test for binary alignment
binalign1 = alntobin(alignment1);
binselection1 = selectalpha(binalign1, [1 3]);
binselection1_exp = alntobin(alnadd(subalign1, subalign3));
utexpect(isequal(binselection1, binselection1_exp), ...
    'selectalpha with binary alignment');

% test with repetitions and order changes
alphaselection3 = selectalpha(alignment1, [2 3 4 2]);
alphaselection3_exp = alnadd(subalign2, alnadd(subalign3, alnadd(subalign4, subalign2)));
utexpect(alncompare(alphaselection3, alphaselection3_exp), ...
    'selectalpha with repetitions and order changes, character alignment');

% test for statistics structure
stats1 = getstats(binalign1);
statsselection1 = selectalpha(stats1, [1 3 2 4 4]);
binselection2 = selectalpha(binalign1, [1 3 2 4 4]);
statsselection1_exp = getstats(binselection2);
utexpect(statscompare(statsselection1, statsselection1_exp), ...
    'selectalpha with statistics structure');

% with freq2
stats1 = getstats(binalign1, 'freq2', true);
statsselection1 = selectalpha(stats1, [1 4 1 3 2]);
binselection2 = selectalpha(binalign1, [1 4 1 3 2]);
statsselection1_exp = getstats(binselection2, 'freq2', true);
utexpect(statscompare(statsselection1, statsselection1_exp), ...
    'selectalpha with statistics structure with ''freq2''');

% using logical mask
binselection3 = selectalpha(binalign1, [false true true false]);
binselection3a = selectalpha(binalign1, [2 3]);
utexpect(isequal(binselection3, binselection3a), ...
    'selectalpha with logical selection');

end