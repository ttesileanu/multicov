function test_getfreq1
% TEST_GETFREQ1 Test getfreq1.m.

nreps = 3;
alphabets = {'protein', 'dna'};
N = 256;
n = 32;
for a = 1:length(alphabets)
    timer = tic;
    alphabet = alphabets{a};
    testres = true;
    letters = alphagetletters(alphabet, 'nogap');
    w = length(letters)*n;
    for i = 1:nreps
        alignment = alngenrandom(N, n, alphabet);
        freq1 = getfreq1(alignment);
        
        freq1exp = zeros(w, 1);
        for j = 1:w
            letter = letters(1 + mod(j-1, length(letters)));
            freq1exp(j) = mean(alignment.data(:, 1 + floor((j-1)/length(letters))) == letter);
        end
        
        if any(abs(freq1 - freq1exp) > eps)
            testres = false;
            break;
        end
    end
    utexpect(testres, ['getfreq1 for ' alphabet], timer);
end

% test with sequence weights and multiple alphabets
alignment1 = alnmake(['AAA' ; 'CCC' ; 'GGC'], 'dna');
alignment2 = alnmake(['AW--W' ; 'DFGYA' ; 'T-SRA'], 'protein');
alignment = alnadd(alignment1, alignment2);
alignment.seqw = [1 ; 1/2 ; 1/4];
freq1 = getfreq1(alignment);
freq1exp = ...
   [1 ; 1/2 ; 1/4 ; 0 ; ...
    1 ; 1/2 ; 1/4 ; 0 ; ...
    1 ; 3/4 ;  0  ; 0 ; ...
%   A    C     D    E   F   G   H   I   K   L   M   N   P   Q   R   S    T    V   W   Y 
    1 ;  0  ; 1/2 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 1/4 ; 0 ; 0 ; 0 ; ...
    0 ;  0  ;  0  ; 0 ;1/2; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ;  0  ; 0 ; 1 ; 0 ; ...
    0 ;  0  ;  0  ; 0 ; 0 ;1/2; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ;1/4;  0  ; 0 ; 0 ; 0 ; ...
    0 ;  0  ;  0  ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ;1/4; 0 ;  0  ; 0 ; 0 ;1/2; ...
   3/4;  0  ;  0  ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ;  0  ; 0 ; 1 ; 0 ...
    ] / sum(alignment.seqw);
utexpect(all(abs(freq1 - freq1exp) < eps), 'getfreq1 with seqw, multiple alphabets');

binalign = alntobin(alignment);
utexpect(all(abs(freq1exp - getfreq1(binalign)) < eps), 'getfreq1 with binary alignment');

% XXX test with binary alignment

% test with stats structure
stats.type = 'stats';
stats.alphabets = {'dna'};
stats.alphawidths = 2;
stats.freq1 = [0.5 ; 0.3 ; 0.1 ; 0.0 ; 0.2 ; 0.2 ; 0.2 ; 0.1];
stats.cmat = stats.freq1*stats.freq1';
stats.refseq.seqdb = 'generic';
stats.refseq.seqid = '';
stats.refseq.map = [1 ; 2];
utexpect(isequal(stats.freq1, getfreq1(stats)), 'getfreq1 with stats structure');

end