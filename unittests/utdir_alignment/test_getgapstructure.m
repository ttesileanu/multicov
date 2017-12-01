function test_getgapstructure
% TEST_GETGAPSTRUCTURE Test getgapstructure.m.

alignment = alngenrandom(256, 64, 'protein');
gaps = (alignment.data == '-');
utexpect(isequal(gaps, getgapstructure(alignment)), 'getgapstructure single alphabet');

pattern.alphabets = {'protein', 'dna', 'protein'};
pattern.alphawidths = [51 ; 3 ; 10];
alignment2 = alngenrandom(512, pattern);
gaps2 = (alignment2.data == '-');
binalign2 = alntobin(alignment2);
utexpect(isequal(gaps2, getgapstructure(alignment2)), ...
    'getgapstructure multi-alphabet');
utexpect(isequal(gaps2, getgapstructure(binalign2)), ...
    'getgapstructure binary');

% test that this does the right thing for gapped alphabets
gapalignment2 = alignment2;
gapalignment2.alphabets{1} = ['gap' alignment2.alphabets{1}];
gapalignment2.alphabets{2} = ['gap' alignment2.alphabets{2}];
gapbinalign2 = alntobin(gapalignment2);
utexpect(isequal(getgapstructure(alignment2), getgapstructure(gapalignment2)) && ...
    isequal(getgapstructure(binalign2), getgapstructure(gapbinalign2)), ...
    'getgapstructure with gapped alphabets');

end