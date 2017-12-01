function test_alnmake
% TEST_ALNMAKE Test alnmake.m.

aln0 = alnmake;
aln0exp.alphabets = {};
aln0exp.alphawidths = [];
aln0exp.data = '';
aln0exp.seqw = [];
aln0exp.refseq = struct('seqdb', {}, 'seqid', {}, 'map', {});
aln0exp.annotations = {};
aln0exp.type = 'character';

utexpect(isequal(aln0, aln0exp), 'alnmake empty alignment');

timer = tic;
letters = '01';
N = 1000;
n = 100;
data = letters(randi(length(letters), N, n));
aln1 = alnmake(data, 'binary');
aln1exp.alphabets = {'binary'};
aln1exp.alphawidths = n;
aln1exp.data = data;
aln1exp.seqw = ones(N, 1);
aln1exp.refseq = struct('seqdb', {'generic'}, 'seqid', {''}, 'map', {(1:n)'});
aln1exp.annotations = repmat({''}, N, 1);
aln1exp.type = 'character';

utexpect(isequal(aln1, aln1exp), 'alnmake binary no annotations', timer);

timer = tic;
letters = '-ACGT';
N = 1200;
n = 64;
data = letters(randi(length(letters), N, n));
annotations = arrayfun(@(i) [int2str(i) '-' int2str(sum(data(i, :) == 'A'))], 1:N, 'uniform', false);
aln2 = alnmake(data, 'dna', annotations);
aln2exp.alphabets = {'dna'};
aln2exp.alphawidths = n;
aln2exp.data = data;
aln2exp.seqw = ones(N, 1);
aln2exp.refseq = struct('seqdb', {'generic'}, 'seqid', {''}, 'map', {(1:n)'});
aln2exp.annotations = annotations;
aln2exp.type = 'character';

utexpect(isequal(aln2, aln2exp), 'alnmake dna with annotations', timer);

end