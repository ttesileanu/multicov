function test_drawfitline
% TEST_DRAWFITLINE Test the non-drawing aspects of drawfitline.m.

x = rand(1000, 1);
y = randn(100, 10);
c = drawfitline(x, y, 'nodraw', true);
utexpect(abs(c - corr(x(:), y(:))) < eps, 'drawfitline correlation coefficient');

cs = drawfitline(x, y, 'nodraw', true, 'corrtype', 'spearman');
utexpect(abs(cs - corr(x(:), y(:), 'type', 'spearman')) < eps, ...
    'drawfitline Spearman correlation');

mocka = 0.13;
mockb = 2.3;
mocky = mocka*x + mockb;
[c, stats] = drawfitline(x, mocky, 'nodraw', true);
tol = 1e-12;
utexpect(abs(1 - c) < tol && abs(stats.a - mocka) < tol && abs(stats.b - mockb) < tol, ...
    'drawfitline coefficients for perfect fit');

[~, stats0] = drawfitline(x, y, 'nodraw', true, 'intercept', 0.5);
utexpect(abs(stats0.b - 0.5) < eps, 'drawfitline fix intercept');

ref1 = [1 ; 2 ; 4 ; 5 ; 6 ; 9 ; 10 ; 11 ; 13 ; 15 ; 17];
ref2 = [2 ; 3 ; 4 ; 5 ; 6 ; 7 ; 8 ; 12 ; 13 ; 14 ; 15 ; 16];
x1 = randn(size(ref1));
x2 = randn(size(ref2));
[c0, stats0] = drawfitline(x1, x2, 'nodraw', true, 'refseq', {ref1, ref2});
subx1 = x1([2 ; 3 ; 4 ; 5 ; 9 ; 10]);
subx2 = x2([1 ; 3 ; 4 ; 5 ; 9 ; 11]);
subc = corr(subx1(:), subx2(:));
subab = polyfit(subx1(:), subx2(:), 1);
utexpect(abs(subc - c0) < eps && all(abs([stats0.a ; stats0.b] - subab(:)) < eps), ...
    'drawfitline with refseq');

end