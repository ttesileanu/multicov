function test_extendfreq
% TEST_EXTENDFREQ Test extendfreq.m.

pattern.alphabets = {'protein'};
pattern.alphawidths = 5;
f0 = rand(20, 1);
freq1 = extendfreq(pattern, f0);
utexpect(all(abs(freq1(:) - repmat(f0, 5, 1)) < eps), ...
    'extendfreq single alphabet, given frequencies');

timer = tic;
pattern.alphabets = {'dna' ; 'protein' ; 'rna'};
pattern.alphawidths = [17 ; 3 ; 10];
[freq1, freq2] = extendfreq(pattern, 'uniform');
freq1exp = [ones(4*17, 1)/5 ; ...
    ones(20*3, 1)/21 ; ...
    ones(4*10, 1)/5];
freq2exp = freq1exp*freq1exp';
binmap = getbinmap(pattern);
for i = 1:length(binmap)
    idxs = binmap{i};
    freq2exp(idxs, idxs) = diag(freq1exp(idxs));
end
utexpect(all(abs(freq1(:) - freq1exp(:)) < eps) && ...
    all(abs(flatten(freq2 - freq2exp)) < eps), ...
    'extendfreq multi alphabet, uniform', timer);

timer = tic;
[gfreq1, gfreq2] = extendfreq(alnincludegaps(pattern), 'uniform');
gfreq1exp = [ones(5*17, 1)/5 ; ...
    ones(21*3, 1)/21 ; ...
    ones(5*10, 1)/5];
gfreq2exp = gfreq1exp*gfreq1exp';
gbinmap = getbinmap(alnincludegaps(pattern));
for i = 1:length(gbinmap)
    idxs = gbinmap{i};
    gfreq2exp(idxs, idxs) = diag(gfreq1exp(idxs));
end
utexpect(all(abs(gfreq1(:) - gfreq1exp(:)) < eps) && ...
    all(abs(flatten(gfreq2 - gfreq2exp)) < eps), ...
    'extendfreq gapped alphabets', timer);

[freq1, freq2] = extendfreq(pattern, {'uniform', 'database', [0.1 0.2 0.3 0.3]});
freq1exp = [ones(4*17, 1)/5 ; ...
    repmat(flatten(getdbbkg('protein')), 3, 1) ; ...
    repmat([0.1 ; 0.2 ; 0.3 ; 0.3], 10, 1)];
freq2exp = freq1exp*freq1exp';
for i = 1:length(binmap)
    idxs = binmap{i};
    freq2exp(idxs, idxs) = diag(freq1exp(idxs));
end
utexpect(all(abs(freq1 - freq1exp) < eps) && ...
    all(abs(flatten(freq2 - freq2exp)) < eps), ...
    'extendfreq multi alphabet, per-alphabet choices', timer);

end