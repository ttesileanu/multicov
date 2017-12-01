function test_selectseq
% TEST_SELECTSEQ Test selectseq.m.

pattern.alphabets = {'dna', 'protein'};
pattern.alphawidths = [33 ; 17];
N = 1024;
alignment = alngenrandom(N, pattern);
alignment.seqw = 0.01 + 0.98*rand(N, 1);
alignment.annotations = arrayfun(@(i) char('a' + randi(26, 1, 10) - 1), ...
    1:N, 'uniform', false);

% test with a permutation
idxs1 = randperm(N);
subalign1 = selectseq(alignment, idxs1);
subalign1_exp = alignment;
subalign1_exp.data = alignment.data(idxs1, :);
subalign1_exp.seqw = alignment.seqw(idxs1);
subalign1_exp.annotations = alignment.annotations(idxs1);

utexpect(isequal(subalign1, subalign1_exp), 'selectseq with permutation, character alignment');

% test with missing sequences & binary alignment
idxs2 = randperm(N);
idxs2 = idxs2(1:floor(N/3));
binalign = alntobin(alignment);
subbinalign2 = selectseq(binalign, idxs2);
subbinalign2_exp = binalign;
subbinalign2_exp.data = binalign.data(idxs2, :);
subbinalign2_exp.seqw = binalign.seqw(idxs2);
subbinalign2_exp.annotations = binalign.annotations(idxs2);

utexpect(isequal(subbinalign2, subbinalign2_exp), ...
    'selectseq with subset of sequences, binary alignment');

% test with repeated indices
idxs3a = randperm(N);
idxs3a = idxs3a(1:floor(N/3));
idxs3b = randperm(N);
idxs3b = idxs3b(1:floor(N/3));
idxs3c = randperm(N);
idxs3c = idxs3c(1:floor(N/3));

idxs3 = [idxs3a(:) ; idxs3b(:) ; idxs3c(:)];
subalign3 = selectseq(alignment, idxs3);

utexpect(isequal(subalign3.data, alignment.data(idxs3, :)) && ...
    isequal(flatten(subalign3.seqw), flatten(alignment.seqw(idxs3))) && ...
    isequal(flatten(subalign3.annotations), flatten(alignment.annotations(idxs3))), ...
    'selectseq with repeated indices');

% test with logical mask
mask = (rand(N, 1) < 0.5);
subalign4 = selectseq(alignment, mask);

utexpect(isequal(subalign4.data, alignment.data(mask, :)) && ...
    isequal(subalign4.seqw, alignment.seqw(mask)) && ...
    isequal(subalign4.annotations, alignment.annotations(mask)), ...
    'selectseq with logical mask');

end