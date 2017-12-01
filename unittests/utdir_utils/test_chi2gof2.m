function test_chi2gof2
% TEST_CHI2GOF2 Test chi2gof2.m.

m1 = [5 7 ; 9 8];
[p1y, chi1y] = chi2gof2(m1, 'yates', true);
[p1n, chi1n] = chi2gof2(m1, 'yates', false);
utexpect(abs(chi1y - 0.049) < 0.001 && ...
    abs(chi1n - 0.358) < 0.001 && ...
    abs(p1y - 0.825) < 0.001 && ...
    abs(p1n - 0.5496) < 0.001, ...
    'chi2gof2 2x2');

m2 = [5 3 12 ; 5 7 15];
[p2y, chi2y] = chi2gof2(m2, 'yates', true);
[p2n, chi2n] = chi2gof2(m2, 'yates', false);
utexpect(all(abs([p2y p2n] - [0.879 0.634]) < 0.001) && ...
    all(abs([chi2y chi2n] - [0.258 0.911]) < 0.001), ...
    'chi2gof2 2x3');

m3 = [10 9 ; 12 3 ; 10 4];
[p3y, chi3y] = chi2gof2(m3, 'yates', true);
[p3n, chi3n] = chi2gof2(m3, 'yates', false);
utexpect(all(abs([p3y p3n] - [0.407 0.220]) < 0.001) && ...
    all(abs([chi3y chi3n] - [1.796 3.027]) < 0.001), ...
    'chi2gof2 3x2');

m4 = [...
     3 15  2 10 13 ; ... 
    12 13 14 10 14 ; ...
     9  6 11 12  4 ; ...
    11 11 13 13  4 ; ...
     3  1  5  6 12 ...
];
[p4y, chi4y] = chi2gof2(m4, 'yates', true);
[p4n, chi4n] = chi2gof2(m4, 'yates', false);
utexpect(all(abs([p4y p4n] - [0.0163 0.0014]) < 0.0001) && ...
    all(abs([chi4y chi4n] - [30.351 38.203]) < 0.001), ...
    'chi2gof2 5x5');

end