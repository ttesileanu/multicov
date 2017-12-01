function test_findminimum
% TEST_FINDMINIMUM Test findminimum.m.

% test with a simple function
fct0 = @(x) norm(x - [1.5 2 pi])^2 + 1;
fct0_der = @(x) (x - [1.5 2 pi]);

timer = tic;
x0 = findminimum(fct0, [0 0 0], 'derivative', fct0_der, 'method', 'descent');
x0_cg = findminimum(fct0, [0 0 0], 'derivative', fct0_der, 'method', 'cg');
utexpect(all(abs(x0 - [1.5 2 pi]) < 0.05), 'findminimum with simple function', timer);

utexpect(abs(x0 - x0_cg) < 0.0001, 'findminimum conjugate-gradient vs. gradient descent');

% test with a more complicated function, but in 1d
fct1 = @(x) sin(4*x) + 0.05*x.^4 + 3*x.^2;
fct1_der = @(x) 4*cos(4*x) + 0.2*x.^3 + 6*x;

timer = tic;
[x1, stats1] = findminimum(fct1, 0, 'derivative', fct1_der, 'rate', 0.1);
utexpect(abs(x1 + 0.2829) < 0.0001, 'findminimum with complicated function, 1d', timer);

utexpect(isstruct(stats1) && all(isfield(stats1, {'result', 'fmin', 'nfeval', 'ndeval'})) && ...
    strcmp(stats1.result, 'converged') && abs(stats1.fmin - fct1(x1)) < 0.0001, ...
    'findminimum convergence statistics');

[~, stats1_mock] = findminimum(fct1, 0, 'derivative', fct1_der, 'maxiter', 1);
utexpect(strcmp(stats1_mock.result, 'failed') && stats1_mock.fmin > stats1.fmin, ...
    'findminimum non-convergence detection, and test of maxiter');

% a final function, this time with many dimensions
fct2 = @(x) norm(x - magic(5))^2;
fct2_der = @(x) 2*(x - magic(5));

timer = tic;
[x2, stats2] = findminimum(fct2, zeros(5), 'derivative', fct2_der, 'tol', 1e-4, 'rate', 0.1);
[x2_better, stats2_better] = findminimum(fct2, zeros(5), 'derivative', fct2_der, 'tol', 1e-5, 'rate', 0.1);
utexpect(all(abs(x2 - magic(5)) < 1e-4), 'findminimum with matrix variable', timer);

utexpect(stats2_better.fmin <= stats2.fmin && ...
    stats2_better.nfeval >= stats2.nfeval && stats2_better.ndeval >= stats2.ndeval && ...
    stats2_better.niter >= stats2.niter && all(abs(x2_better(:) - x2(:)) < 1e-4), ...
    'findminimum tolerance test');

% XXX need more complex tests, these finish in a 2-3 iterations...

end