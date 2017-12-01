function test_getfreq2
% TEST_GETFREQ2 Test getfreq2.m.

nreps = 2;
alphabets = {'protein', 'dna'};
N = 128;
n = 7;
for a = 1:length(alphabets)
    timer = tic;
    alphabet = alphabets{a};
    testres = true;
    letters = alphagetletters(alphabet, 'nogap');
    w = length(letters)*n;
    for i = 1:nreps
        alignment = alngenrandom(N, n, alphabet);
        freq2 = getfreq2(alignment);
        
        freq2exp = zeros(w, w);
        for p = 1:w
            letter1 = letters(1 + mod(p-1, length(letters)));
            freq2exp(p, p) = mean(alignment.data(:, 1 + floor((p-1)/length(letters))) == letter1);
            for q = p+1:w
                letter2 = letters(1 + mod(q-1, length(letters)));
                freq2exp(p, q) = mean(...
                    alignment.data(:, 1 + floor((p-1)/length(letters))) == letter1 & ...
                    alignment.data(:, 1 + floor((q-1)/length(letters))) == letter2);
                freq2exp(q, p) = freq2exp(p, q);
            end
        end
        
        if any(abs(freq2(:) - freq2exp(:)) > eps)
            testres = false;
            break;
        end
    end
    utexpect(testres, ['getfreq2 for ' alphabet], timer);
end

% test with sequence weights and multiple alphabets
pattern.alphabets = {'dna', 'protein', 'rna'};
pattern.alphawidths = [5 ; 3 ; 4];
alignment2 = alngenrandom(N, pattern);
alignment2.seqw = rand(N, 1);
w = sum(arrayfun(...
    @(i) length(alphagetletters(pattern.alphabets{i}, 'nogap'))*pattern.alphawidths(i), ...
    1:length(pattern.alphabets)));
freq2exp = zeros(w, w);
binmap = getbinmap(pattern);
cumw = cumsum(pattern.alphawidths);
for i = 1:sum(pattern.alphawidths)
    alpha1 = pattern.alphabets{find(i <= cumw, 1)};
    letters1 = alphagetletters(alpha1, 'nogap');
    for ci = 1:length(letters1)
        lett1 = letters1(ci);
        mask1 = (alignment2.data(:, i) == lett1);
        for j = 1:sum(pattern.alphawidths)
            alpha2 = pattern.alphabets{find(j <= cumw, 1)};
            letters2 = alphagetletters(alpha2, 'nogap');
            for cj = 1:length(letters2)
                lett2 = letters2(cj);
                mask2 = (alignment2.data(:, j) == lett2);
                freq2exp(binmap{i}(ci), binmap{j}(cj)) = sum((mask1 & mask2) .* alignment2.seqw);
            end
        end
    end
end
freq2exp = freq2exp / sum(alignment2.seqw);
freq2 = getfreq2(alignment2);
utexpect(all(abs(freq2(:) - freq2exp(:)) < eps), 'getfreq2 with seqw, multiple alphabets');

binalign = alntobin(alignment2);
freq2b = getfreq2(binalign);
utexpect(all(abs(freq2(:) - freq2b(:)) < eps), 'getfreq2 with binary alignment');

% tests with stats structure
stats.type = 'stats';
stats.alphabets = {'dna'};
stats.alphawidths = 2;
stats.freq1 = [0.5 ; 0.3 ; 0.1 ; 0.0 ; 0.2 ; 0.2 ; 0.2 ; 0.1];
stats.freq2 = stats.freq1*stats.freq1';
stats.cmat = zeros(size(stats.freq2));
stats.refseq.seqdb = 'generic';
stats.refseq.seqid = '';
stats.refseq.map = [1 ; 2];
utexpect(isequal(stats.freq2, getfreq2(stats)), 'getfreq2 with stats structure with freq2 field');

stats = rmfield(stats, 'freq2');
utexpect(isequal(stats.freq1*stats.freq1', getfreq2(stats)), 'getfreq2 with stats structure without freq2 field');

end