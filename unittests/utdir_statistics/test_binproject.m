function test_binproject
% TEST_BINPROJECT Test binproject.m.

pattern.alphabets = {'dna', 'rna', 'protein'};
pattern.alphawidths = [32 ; 10 ; 7];
alignment = alngenrandom(256, pattern);
binalign = alntobin(alignment);
vec = 0.01 + 0.99*rand(size(binalign.data, 2), 1);
binproj = binproject(binalign, vec);
projexp = zeros(size(alignment.data));
binmap = getbinmap(alignment);
for i = 1:length(binmap)
    subvec = vec(binmap{i});
    subvec = subvec / norm(subvec);
    projexp(:, i) = binalign.data(:, binmap{i})*subvec(:);
end
utexpect(all(abs(projexp(:) - binproj.data(:)) < eps), ...
    'binproject');

end