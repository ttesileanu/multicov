function test_comparesectors
% TEST_COMPARESECTORS Test comparesectors.m.

% a simple test
k = 5;
n = 100;
s = 20;
secs1 = cell(1, k);
for i = 1:k
    secs1{i} = randperm(n, s);
end
rndperm = randperm(k);
secs2 = secs1(rndperm);
[ssim, smap] = comparesectors(secs2, secs1);
utexpect(all(abs(1 - ssim) < eps) && all(rndperm(:) == smap(:)), ...
    'comparesectors permutation test');

% test with strings
secsb_1 = {'5', '12', '15', '7', '8'};
secsb_2 = {{'5', '11', '15', '7'}, [1, 2, 3, 4], {'5', '13', '14', '9', '10'}};
[ssim, smap, smat] = comparesectors(secsb_1, secsb_2, 'method', 'first');
smat_exp = [3/5 0 1/5];
smap_exp = 1;
ssim_exp = smat_exp(smap_exp);
utexpect(all(abs(ssim_exp - ssim) < eps) && ...
    all(abs(smap_exp - smap) < eps) && ...
    all(abs(smat_exp - smat) < eps), ...
    'comparesectors various shapes, different method');

end