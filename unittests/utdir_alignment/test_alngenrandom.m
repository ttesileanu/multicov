function test_alngenrandom
% TEST_ALNGENRANDOM Test aln_genrandom.m.

N = 256;
n = 32;
nreps = 8;

pattern1.alphabets = {'protein'};
pattern1.alphawidths = n;

% to ensure reproducibility, always seed the random number generator in the
% same way
timer = tic;
rng('default');
testres = agr_singletest(ones(20*n, 1)/21, nreps, N, pattern1);
utexpect(testres, 'alngenrandom uniform, from pattern', timer);

timer = tic;
rng('default');
testres = agr_singletest(ones(4*n, 1)/5, nreps, N, n, 'dna');
utexpect(testres, 'alngenrandom uniform, from alphabet', timer);

timer = tic;
rng('default');
sf_fct0 = @(i) sin(i)^2*exp(linspace(0.01, 1, 21));
normalize_fct = @(v) v/sum(v);
chop_gap = @(v) v(2:end);
sf_fct = @(i) chop_gap(normalize_fct(sf_fct0(i)));
freq1 = cell2mat(arrayfun(sf_fct, 1:n, 'uniform', false));
freq1 = freq1(:);
alignment1 = alngenrandom(N, freq1, pattern1);
testres = agr_singletest(freq1, nreps, N, freq1, pattern1);
utexpect(testres, 'alngenrandom non-uniform, from pattern', timer);

timer = tic;
utexpect(~agr_singletest(ones(20*n, 1)/21, nreps, N, freq1, pattern1), ...
    'alngenrandom non-uniform detection', timer);

timer = tic;
rng('default');
alignment2 = alngenrandom(N, freq1, 'protein');
testres = agr_singletest(freq1, nreps, N, freq1, 'protein');
utexpect(testres && alncompare(alignment1, alignment2), 'alngenrandom non-uniform, from alphabet', timer);

exp_refseq.seqdb = 'generic';
exp_refseq.seqid = '';
exp_refseq.map = (1:n)';
utexpect(isequal(alignment1.refseq, exp_refseq), 'alngenrandom refseq');

pattern2.alphabets = {'protein', 'dna'};
pattern2.alphawidths = [n 5];

timer = tic;
rng('default');
testres = agr_singletest([ones(20*n, 1)/21 ; ones(4*5, 1)/5], nreps, N, pattern2);
utexpect(testres, 'alngenrandom multi-alphabet', timer);

% test that this does the right thing for gapped alphabets
gappattern3.alphabets = {'gapprotein', 'dna', 'gaprna'};
gappattern3.alphawidths = [15 ; 23 ; 8];
gappattern3.refseq = struct('seqdb', {'uniprot', 'generic', 'trnace'}, ...
    'seqid', {'fake123', 'none', 'fake234'}, ...
    'map', {5:2:33, 2:24, 5:3:26});
alignment3 = alngenrandom(N, gappattern3);
letts_gprot = alphagetletters('gapprotein');
letts_dna = alphagetletters('dna');
letts_grna = alphagetletters('gaprna');
utexpect(all(flatten(alignment3.data(:, 1:15) ~= letts_gprot(1))) && ...
    all(flatten(alignment3.data(:, 38:end) ~= letts_grna(1))) && ...
    any(flatten(alignment3.data(:, 16:38) == letts_dna(1))) && ...
    isequal(alignment3.refseq, gappattern3.refseq), ...
    'alngenrandom gapped alphabets, with pattern (including refseq)');

alignment4 = alngenrandom(N, n, 'gapmulti5');
letts_gmulti5 = alphagetletters('gapmulti5');
utexpect(all(alignment4.data(:) ~= letts_gmulti5(1)), ...
    'alngenrandom gapped single alphabet');

end