function test_linesearch
% TEST_LINESEARCH Test linesearch.m.

timer = tic;
fct1 = @(x) sin(x).^2 + 0.5*cos(0.6*x - 0.3);
[x1, stats] = linesearch(fct1, 2, 3, 4);
utexpect(abs(x1 - 3.2931) < 2e-4, 'linesearch with known brackets', timer);

utexpect(isstruct(stats) && all(isfield(stats, {'result', 'fmin', 'niter', 'neval'})) && ...
    strcmp(stats.result, 'converged') && abs(stats.fmin - fct1(x1)) < 1e-4, ...
    'linesearch convergence statistics check');

[~, stats_short] = linesearch(fct1, 2, 3, 4, 'maxiter', 3);
utexpect(strcmp(stats_short.result, 'failed'), 'linesearch detection of failure to converge');

timer = tic;
[x2, stats2] = linesearch(fct1, 2.1, 2.8, 'tol', 1e-6);
utexpect(abs(x2 - x1) < 2e-3, 'linesearch with known direction', timer);

timer = tic;
x2b = linesearch(fct1, 2.8, 2.1, 'tol', 1e-6);
utexpect(abs(x2b - x2) < 2e-3, 'linesearch with known direction, reverse sign', timer);

timer = tic;
[~, stats2_better] = linesearch(fct1, 2.1, 2.8, 'tol', 1e-8);
utexpect(stats2_better.fmin <= stats2.fmin && stats2_better.niter > stats2.niter && ...
    stats2_better.neval > stats2.neval, 'linesearch changing the tolerance', timer);

timer = tic;
fct2 = @(x) sin(4*x) + 0.05*x.^4 + 3*x.^2;

x2 = linesearch(fct2, 0, 0.01);
utexpect(abs(x2 + 0.2829) < 0.0001, 'linesearch with complicated function', timer);

end