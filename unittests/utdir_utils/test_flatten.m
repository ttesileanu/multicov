function test_flatten
% TEST_FLATTEN Test flatten.m

n = 32;
v1 = randn(n, 1);
v2 = randn(1, n);
utexpect(all(abs(flatten(v1) - v1) < eps) && all(abs(flatten(v2) - v2') < eps), ...
    'flatten for vectors');

m = randn(n, n);
utexpect(all(abs(flatten(m) - m(:)) < eps), 'flatten for matrices');

t = randn(5, 7, 13, 3, 4);
utexpect(all(abs(flatten(t) - t(:)) < eps), 'flatten for higher-dimensional array');

end