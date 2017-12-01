function test_eigsorted
% TEST_EIGSORTED Test eigsorted.m.

% XXX Many of these calculations give errors much larger than eps... what
% exactly are eig's guarantees?
tol = 1e-12;

% simple check for diagonal matrix
l = randn(16, 1);
m = diag(l);
utexpect(all(abs(eigsorted(m) - sort(l, 'descend')) < tol), 'eigsorted eigenvals diagonal matrix');

% check that the ordering works well
m2 = randn(32, 32);
m2 = m2 + m2';
[v2, d2] = eigsorted(m2);
[v2exp, d2exp] = eig(m2);
d2exp = diag(d2exp);
[d2exp, ord2] = sort(d2exp, 'descend');
v2exp = v2exp(:, ord2);
utexpect(all(abs(d2(:) - d2exp(:)) < tol) && ...
    all(arrayfun(@(i) (all(abs(v2(:, i) - v2exp(:, i)) < tol) || ...
    all(abs(v2(:, i) + v2exp(:, i)) < tol)), 1:size(m2, 1))), ...
    'eigsorted evals&evecs for random matrix');

% check that we're assigning the signs of the eigenvectors as promised
testres = true;
for i = 1:size(m2, 1)
    onev = v2(:, i);
    [~, mx] = max(abs(onev));
    if onev(mx) < 0
        testres = false;
        break;
    end
end
utexpect(testres, 'eigsorted eigenvector signs');

% restricting to top eigenvectors/eigenvalues
k = 5;
dchop2a = eigsorted(m2, k);
utexpect(all(abs(dchop2a - d2(1:k)) < tol), ...
    'eigsorted top k eigenvals');

[vchop2, dchop2] = eigsorted(m2, k);
utexpect(all(abs(dchop2a - dchop2) < tol) && all(all(abs(vchop2 - v2(:, 1:k)) < tol)), ...
    'eigsorted top k evals&evecs');

end