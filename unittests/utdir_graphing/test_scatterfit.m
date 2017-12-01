function test_scatterfit
% TEST_SCATTERFIT Test that scatterfit.m returns the same results as
% drawfitline.m.

ntests = 5;
testres = true;
tol = 1e-12;
for i = 1:ntests
    x = randn(1000, 1);
    y = rand(100, 10);
    [c1, stats1] = drawfitline(x, y, 'nodraw', true);
    [c2, stats2] = scatterfit(x, y, 'nodraw', true);
    testres = testres & (...
        abs(c1 - c2) < tol && ...
        all(abs([stats1.a stats1.b] - [stats2.a stats2.b]) < tol) && ...
        all(all(abs(stats1.ci - stats2.ci) < tol)) && ...
        all(abs(stats1.residuals - stats2.residuals) < tol));
end
utexpect(testres, 'scatterfit matches drawfitline');

end