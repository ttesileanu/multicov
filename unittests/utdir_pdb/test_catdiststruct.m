function test_catdiststruct
% TEST_CATDISTSTRUCT Test catdiststruct.m.

symmetrize = @(m) (m + m')/2;
n1 = 16;
n2 = 24;

dist1.distmat = symmetrize(randn(n1, n1));
dist1.validmask = (symmetrize(rand(n1, n1)) < 0.9);
dist1.refseq.seqdb = 'a';
dist1.refseq.id = '1';
dist1.refseq.map = 1:2:(2*n1-1);

dist2.distmat = symmetrize(randn(n2, n2));
dist2.validmask = (symmetrize(rand(n2, n2)) < 0.85);
dist2.refseq.seqdb = 'b';
dist2.refseq.id = '2';
dist2.refseq.map = 1:n2;

dist = catdiststruct(dist1, dist2);
utexpect(isequal(dist.distmat(1:n1, 1:n1), dist1.distmat) && ...
    isequal(dist.distmat((n1+1):end, (n1+1):end), dist2.distmat), ...
    'catdiststruct within-block distances correct');

utexpect(all(abs(flatten(dist.distmat(1:n1, (n1+1):end))) < eps) && ...
    all(abs(flatten(dist.distmat - dist.distmat')) < eps) && ...
    ~any(flatten(dist.validmask(1:n1, (n1+1):end))) && ...
    all(flatten(dist.validmask == dist.validmask')), ...
    'catdiststruct distmat symmetric and distances between blocks set to invalid zero');

utexpect(length(dist.refseq) == 2 && isequal(dist.refseq(1), dist1.refseq) && ...
    isequal(dist.refseq(2), dist2.refseq), 'catdiststruct correct refseq');

dist_again = catdiststruct(dist, dist1);
dist_exp = zeros(2*n1 + n2, 2*n1 + n2);
dist_exp(1:n1, 1:n1) = dist1.distmat;
dist_exp((n1+1):(n1+n2), (n1+1):(n1+n2)) = dist2.distmat;
dist_exp((n1+n2+1):end, (n1+n2+1):end) = dist1.distmat;
utexpect(all(abs(dist_again.distmat(:) - dist_exp(:)) < eps), ...
    'catdiststruct two-element plus one-element structure correct distances');

valid_exp = [...
    dist1.validmask false(n1, n1 + n2) ; ...
    false(n2, n1) dist2.validmask false(n2, n1) ; ...
    false(n1, n1) false(n1, n2) dist1.validmask];
utexpect(all(valid_exp(:) == dist_again.validmask(:)), ...
    'catdiststruct two-element plus one-element structure correct valid mask');

utexpect(length(dist_again.refseq) == 3 && ...
    isequal(dist_again.refseq(1), dist1.refseq) && ...
    isequal(dist_again.refseq(2), dist2.refseq) && ...
    isequal(dist_again.refseq(3), dist1.refseq), ...
    'catdiststruct two-element plus one-element structure correct refseq');

end