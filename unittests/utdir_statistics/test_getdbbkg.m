function test_getdbbkg
% TEST_GETDBBKG Test getdbbkg.m.

alphabets = alphagetletters('list');
tol = 1e-12;
for i = 1:length(alphabets)
    bkgfreq = getdbbkg(alphabets{i});
    nletters = length(alphagetletters(alphabets{i}, 'nogap'));
    utexpect(all(bkgfreq >= max(0, 1/(5*nletters)) & bkgfreq <= min(1, 5/nletters)) && ...
        abs(sum(bkgfreq) - 1) < tol, ['getdbbkg ' alphabets{i}]);
    
    gapalphabet = ['gap' alphabets{i}];
    bkgfreq1 = getdbbkg(gapalphabet);
    setdbbkg(gapalphabet, randn(size(bkgfreq1)));
    bkgfreq1a = getdbbkg(gapalphabet);
    bkgfreq0a = getdbbkg(alphabets{i});
    setdbbkg(gapalphabet, []);
    
    utexpect(all(abs(bkgfreq0a(:) - bkgfreq(:)) < eps) && all(abs(bkgfreq1(:)' - [0 bkgfreq(:)']) < eps) && ...
        any(abs(bkgfreq1(:) - bkgfreq1a(:)) > eps), ['getdbbkg ' gapalphabet]);
end

end