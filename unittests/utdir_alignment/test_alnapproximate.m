function test_alnapproximate
% TEST_ALNAPPROXIMATE Test alnapproximate.m.

% test binary approximation
N = 256;
n = 32;
alignment0 = alngenrandom(N, n, 'protein');
[alnbinapprox0, binmaps0] = alnapproximate(alignment0, 2);
cons0 = getconsensus(alignment0);
considxs0 = getconsensus(alignment0, 'indices', true);
binalign0 = alntobin(alignment0);

utexpect(strcmp(cons0, binmaps0) && isequal(alnbinapprox0.data == '1', ...
    binalign0.data(:, considxs0)), ...
    'alnapproximate single alphabet, binary approx, including ''maps'' test');

% test non-binary approximation
k = 7;
[alnkapprox0, kmaps0] = alnapproximate(alignment0, k);
freq1_0 = getfreq1(alignment0);
letters0 = alphagetletters(alignment0.alphabets{1}, 'nogap');
multiletters0 = alphagetletters(['multi' int2str(k)]);
kmaps0_exp = repmat(blanks(n), k-1, 1);
kapproxdata0_exp = repmat('0', N, n);
binmap0 = getbinmap(alignment0);
for i = 1:alignment0.alphawidths
    [~, sortorder] = sort(freq1_0(binmap0{i}), 'descend');
    kmaps0_exp(:, i) = letters0(sortorder(1:k-1));
    for j = 1:k-1
        mask = (alignment0.data(:, i) == kmaps0_exp(j, i));
        kapproxdata0_exp(mask, i) = multiletters0(j+1);
    end
end

utexpect(isequal(alnkapprox0.alphabets, {['multi' int2str(k)]}) && ...
    isequal(alnkapprox0.data, kapproxdata0_exp) && isequal(kmaps0, kmaps0_exp), ...
    ['alnaproximate single alphabet, ' int2str(k) '-ary approximation']);

% test multiple alphabets & sequence weights
pattern1.alphabets = {'rna', 'multi6', 'protein'};
pattern1.alphawidths = [7 ; 3 ; 4];
alignment1 = alngenrandom(N, pattern1);
alignment1.seqw = 0.01 + 0.98*rand(N, 1);

p = 5;
alnpapprox1 = alnapproximate(alignment1, p);
papproxdata1_exp = repmat('0', N, sum(pattern1.alphawidths));
multiletters1 = alphagetletters(['multi' int2str(p)]);
freq1_1 = getfreq1(alignment1);
binmap1 = getbinmap(alignment1);
i = 1;
for a = 1:length(alignment1.alphabets)
    letters1 = alphagetletters(alignment1.alphabets{a}, 'nogap');
    for ai = 1:alignment1.alphawidths(a)
        [~, sortorder] = sort(freq1_1(binmap1{i}), 'descend');
        crtmap = letters1(sortorder(1:p-1));
        for j = 1:p-1
            mask = (alignment1.data(:, i) == crtmap(j));
            papproxdata1_exp(mask, i) = multiletters1(j+1);
        end
        i = i + 1;
    end
end

utexpect(isequal(alignment1.alphawidths, alnpapprox1.alphawidths) && ...
    isequal(alnpapprox1.alphabets(:), repmat({['multi' int2str(p)']}, length(alignment1.alphabets), 1)) && ...
    isequal(alnpapprox1.seqw(:), alignment1.seqw(:)) && ...
    isequal(alnpapprox1.annotations(:), alignment1.annotations(:)) && ...
    isequal(alnpapprox1.refseq(:), alignment1.refseq(:)) && ...
    isequal(alnpapprox1.data, papproxdata1_exp), ...
    'alnapproximate multiple alphabest, seqw, other fields');

% test treating gaps like any amino acid
alnpapprox1g = alnapproximate(alignment1, p, 'gapszero', false);
papproxdata1g_exp = repmat('0', N, sum(pattern1.alphawidths));
freq1g_1 = getfreq1(includegaps(alignment1));
binmap1g = getbinmap(includegaps(alignment1));
i = 1;
for a = 1:length(alignment1.alphabets)
    letters1g = alphagetletters(alignment1.alphabets{a});
    for ai = 1:alignment1.alphawidths(a)
        [~, sortorder] = sort(freq1g_1(binmap1g{i}), 'descend');
        crtmap = letters1g(sortorder(1:p-1));
        for j = 1:p-1
            mask = (alignment1.data(:, i) == crtmap(j));
            papproxdata1g_exp(mask, i) = multiletters1(j+1);
        end
        i = i + 1;
    end
end

utexpect(isequal(alnpapprox1g.data, papproxdata1g_exp), ...
    'alnapproximate with ''gapszero'' false');

% test handling of alphabets with fewer than n characters
had_exception = false;
try
    alnapproximate(alignment1, 7);
catch me
    had_exception = true;
end

exp_eprefix = 'alnapproximate:';
utexpect(had_exception && length(me.identifier) > length(exp_eprefix) && ...
    strcmp(me.identifier(1:length(exp_eprefix)), exp_eprefix), ...
    'alnapproximate error when n > number of characters in some alphabet');

end