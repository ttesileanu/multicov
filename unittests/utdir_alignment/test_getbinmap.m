function test_getbinmap
% TEST_GETBINMAP Test getbinmap.m.

pattern.alphabets = {'protein' ; 'dna'};
pattern.alphawidths = [15 ; 32];
map = getbinmap(pattern);

n = sum(pattern.alphawidths);
expmap = cell(n, 1);
crt = 1;
crtbin = 1;
for a = 1:length(pattern.alphabets)
    w = pattern.alphawidths(a);
    nletts = length(alphagetletters(pattern.alphabets{a}, 'nogap'));
    for i = 1:w
        nextbin = crtbin + nletts;
        expmap{crt} = (crtbin:nextbin-1);
        crtbin = nextbin;
        crt = crt + 1;
    end
end
utexpect(isequal(map, expmap), 'getbinmap');

end