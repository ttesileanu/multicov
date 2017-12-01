function test_getsectors
% TEST_GETSECTORS Test getsectors.m.

n = 128;
k = 6;
comps1 = randn(n, k);
secs1_1 = getsectors(comps1, 'cutoff', 1, 'tails', 2, 'method', 'magnitude');
secs1_1_exp = arrayfun(@(i) find(abs(comps1(:, i)) >= 1), 1:k, 'uniform', false);
compare_secs = @(s1, s2) numel(s1) == numel(s2) && all(cellfun(@(a, b) isequal(sort(a(:)), sort(b(:))), s1(:), s2(:)));
utexpect(compare_secs(secs1_1.secs, secs1_1_exp) && isequal(comps1, secs1_1.comps), ...
    'getsectors from comps, by magnitude, two-tails');

sec_map = 3:2:(1+2*128);
secs1_2 = getsectors(comps1, 'cutoff', 1, 'method', 'range', 'tails', 1, ...
    'fixmu', 0, 'refseq', struct('map', sec_map));
sigmas = arrayfun(@(i) norm(comps1(:, i))/sqrt(n-1), 1:k);
%sigmas = std(comps1, 0, 1);
secs1_2_exp = arrayfun(@(i) find(comps1(:, i) >= sigmas(i)), ...
    1:k, 'uniform', false);
secs1_2_mapped_exp = {cellfun(@(v) sec_map(v), secs1_2_exp, 'uniform', false)};
utexpect(compare_secs(secs1_2.secs, secs1_2_exp) && ...
    compare_secs(secs1_2.secs_mapped{1}, secs1_2_mapped_exp{1}), ...
    'getsectors from comps, by range, one-tail, with refseq');

timer = tic;
pattern.alphabets = {'rna', 'protein'};
pattern.alphawidths = [12 ; 7];
alignment = alngenrandom(1024, pattern);
alignment.refseq(1).map = arrayfun(@int2str, 1:12, 'uniform', false);
alignment.refseq(2).map = arrayfun(@int2str, 1:7, 'uniform', false);
sca = getsca(alignment);
secs2_1 = getsectors({sca, 3}, 'method', {'cdf', 'pdf', 'count'}, ...
    'distrib', {'normal', 'normal', 'none'}, 'cutoff', [0.9, 0.1, 6], ...
    'tails', [1, 2, 1]);
[evecs, ~] = eigsorted(sca.cmat);
fit1 = fitdist(evecs(:, 1), 'normal');
fit2 = fitdist(evecs(:, 2), 'normal');
cdf1x = linspace(min(evecs(:, 1)), max(evecs(:, 1)), 4096);
cdf1y = cdf(fit1, cdf1x);
thresh1 = cdf1x(find(cdf1y >= 0.9, 1));
pdf2x = linspace(min(evecs(:, 2)), max(evecs(:, 2)), 4096);
pdf2y = pdf(fit2, pdf2x);
thresh2a = pdf2x(find(pdf2y > 0.1, 1, 'first'));
thresh2b = pdf2x(find(pdf2y <= 0.1, 1, 'last'));
[~, ord3] = sort(evecs(:, 3), 'descend');
secs2_1_exp = {...
    find(evecs(:, 1) >= thresh1), ...
    find(evecs(:, 2) <= thresh2a | evecs(:, 2) >= thresh2b), ...
    ord3(1:6) ...
};
secs2_1_mapped_exp = arrayfun(@(i) cell(1, 3), 1:2, 'uniform', false);
for i = 1:3
    sec = secs2_1_exp{i};
    sec_rna = sec(sec <= 12);
    sec_prot = sec(sec > 12) - 12;
    secs2_1_mapped_exp{1}{i} = alignment.refseq(1).map(sec_rna);
    secs2_1_mapped_exp{2}{i} = alignment.refseq(2).map(sec_prot);
end
utexpect(compare_secs(secs2_1.secs, secs2_1_exp) && ...
    compare_secs(secs2_1.secs_mapped{1}, secs2_1_mapped_exp{1}) && ...
    compare_secs(secs2_1.secs_mapped{2}, secs2_1_mapped_exp{2}), ...
    'getsectors from SCA, with per-sector methods', timer);

stds = [0.2 0.1 0.25];
secs2_2 = getsectors({sca, 3}, 'overlaps', false, 'method', 'range', ...
    'fixmu', 0, 'fixstd', stds, 'cutoff', 1, 'tails', 2);
masks = bsxfun(@ge, abs(evecs(:, 1:3)), stds);
secs2_2_exp = {...
    find(masks(:, 1) & ~masks(:, 2) & ~masks(:, 3)), ...
    find(~masks(:, 1) & masks(:, 2) & ~masks(:, 3)), ...
    find(~masks(:, 1) & ~masks(:, 2) & masks(:, 3)) ...
};
utexpect(compare_secs(secs2_2.secs, secs2_2_exp), ...
    'getsectors from SCA, overlaps disallowed');

evec_choice = [2 4 7];
compsx = evecs(:, evec_choice);
secs3_1 = getsectors(compsx, 'method', 'count', 'cutoff', 6, 'refseq', alignment.refseq);
secs3_2 = getsectors({sca, evec_choice}, 'method', 'count', 'cutoff', 6);
utexpect(compare_secs(secs3_1.secs, secs3_2.secs) && ...
    compare_secs(secs3_1.secs_mapped{1}, secs3_2.secs_mapped{1}) && ...
    compare_secs(secs3_1.secs_mapped{2}, secs3_2.secs_mapped{2}), ...
    'getsectors from SCA, arbitrary selection of evecs');

end