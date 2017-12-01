function test_setcompare
% TEST_SETCOMPARE Test setcompare.m.

n = 25;
k = 10;
s1 = randperm(n, k);
s2 = randperm(n, k);
nint = length(intersect(s1, s2));
nunion = length(union(s1, s2));
m = min(length(s1), length(s2));
M = max(length(s1), length(s2));

du = nint/nunion;
utexpect(abs(du - setcompare(s1, s2, 'method', 'union')) < eps, ...
    'setcompare method ''union''');

dm = nint/m;
utexpect(abs(dm - setcompare(s1, s2, 'method', 'min')) < eps, ...
    'setcompare method ''min''');

dM = nint/M;
utexpect(abs(dM - setcompare(s1, s2, 'method', 'max')) < eps, ...
    'setcompare method ''max''');

d1 = nint/length(s1);
utexpect(abs(d1 - setcompare(s1, s2, 'method', 'first')) < eps, ...
    'setcompare method ''first''');

d2 = nint/length(s2);
utexpect(abs(d2 - setcompare(s1, s2, 'method', 'second')) < eps, ...
    'setcompare method ''second''');

da = 2*nint/(m + M);
utexpect(abs(da - setcompare(s1, s2, 'method', 'mean')) < eps, ...
    'setcompare method ''mean''');

end