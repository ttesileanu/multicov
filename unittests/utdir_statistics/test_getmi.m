function test_getmi
% TEST_GETMI Test getmi.m.

timer = tic;
pattern.alphabets = {'binary', 'multi6'};
pattern.alphawidths = [8 ; 11];
N = 256;
alignment = alngenrandom(N, pattern);
binalign = alntobin(alignment);
stats = getstats(alignment);
mimat1 = getmi(alignment);

gapstats = statsincludegaps(stats);
binmap = getbinmap(gapstats);
n = length(binmap);
mimat_exp = zeros(n);
freq2 = gapstats.cmat + gapstats.freq1(:)*gapstats.freq1(:)';
for i = 1:n
    idxsi = binmap{i};
    freqi = gapstats.freq1(idxsi);
    for j = 1:n
        idxsj = binmap{j};
        freqj = gapstats.freq1(idxsj);
        block = freq2(idxsi, idxsj);
        
        fprod = freqi(:)*freqj(:)';
        
        mask = (fprod > eps & block > eps);
        
        mimat_exp(i, j) = sum(sum(block(mask).*log(block(mask)./fprod(mask))));
    end
end
utexpect(all(abs(mimat1(:) - mimat_exp(:)) < eps), ...
    'getmi from character alignment', timer);

timer = tic;
mimat1a = getmi(binalign);
utexpect(all(abs(mimat1a(:) - mimat1(:)) < eps), ...
    'getmi from binary alignment', timer);

timer = tic;
mimat1b = getmi(stats);
utexpect(all(abs(mimat1b(:) - mimat1a(:)) < eps), ...
    'getmi from statistics structure', timer);

timer = tic;
gapalignment = alnincludegaps(alignment);
gapbinalign = alnincludegaps(binalign);
gapmimat1 = getmi(gapalignment);
gapmimat1a = getmi(gapbinalign);
gapmimat1b = getmi(gapstats);
utexpect(all(abs(gapmimat1(:) - mimat1(:)) < eps) && ...
    all(abs(gapmimat1a(:) - mimat1a(:)) < eps) && ...
    all(abs(gapmimat1b(:) - mimat1b(:)) < eps), ...
    'getmi from gapped structures', timer);

timer = tic;
mockstruct.freq1 = gapstats.freq1;
mockstruct.freq2 = freq2;
mockstruct.alphabets = gapstats.alphabets;
mockstruct.alphawidths = gapstats.alphawidths;
mimock = getmi(mockstruct);
utexpect(all(abs(mimock(:) - mimat1(:)) < eps), ...
    'getmi from gapped stats-like structure', timer);

timer = tic;
mockstruct = rmfield(mockstruct, 'freq2');
mockstruct.cmat = freq2 - mockstruct.freq1(:)*mockstruct.freq1(:)';
mimock = getmi(mockstruct);
utexpect(all(abs(mimock(:) - mimat1(:)) < eps), ...
    'getmi from gapped stats-like structure with cmat but no freq2', timer);

timer = tic;
mockstruct2.freq1 = stats.freq1;
mockstruct2.cmat = stats.cmat;
mockstruct2.alphabets = stats.alphabets;
mockstruct2.alphawidths = stats.alphawidths;
mimock2 = getmi(mockstruct2);
utexpect(all(abs(mimock2(:) - mimat1(:)) < eps), ...
    'getmi from ungapped stats-like structure with cmat but no freq2', timer);

end