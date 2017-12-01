function test_getnth
% TEST_GETNTH Test getnth.m.

d = getnth(2, @eig, [0 1 ; -1 0]);
utexpect(norm(d - diag(diag(d))) < eps && ...
    norm(diag(d) - [1i ; -1i]) < eps, 'getnth function');

m = [1 2 ; 3 4];
utexpect(getnth(2, m(2, :)) == 4, 'getnth vector');

end