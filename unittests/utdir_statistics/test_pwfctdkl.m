function test_pwfctdkl
% TEST_PWFCTDKL Test pwfctdkl.m.

% don't want to test the exact background distribution, but let's test that
% whatever it is, it is reasonable, and is used consistently across sites
pattern1.alphabets = {'protein'};
pattern1.alphawidths = 16;
% avoid 0s and 1s
% make sure we include a probability for gaps, but not store it in freq1
normvec = @(v) v(2:end)/sum(v);
freq1 = cell2mat(arrayfun(@(i) normvec(0.01 + 0.98*rand(21, 1)), (1:pattern1.alphawidths)', 'uniform', false));

posw = pwfctdkl(freq1, pattern1);

% invert to find out what background distribution was used
tmp = exp(posw - log(freq1./(1 - freq1))); % == (1 - bkgfreq1) ./ bkgfreq1
bkgfreq1 = 1./(1+tmp);
bkgfreq1 = bkgfreq1(:);

% check that all sites have the same background distribution
utexpect(all(abs(bkgfreq1 - repmat(bkgfreq1(1:20), pattern1.alphawidths, 1)) < eps), ...
    'pwfctdkl single alphabet, site consistency');

% check that the distribution makes sense
utexpect(all(bkgfreq1 > 0.01 & bkgfreq1 < 0.1), 'pwfctdkl protein, reasonable bkgfreq');

% multiple alphabets
pattern2.alphabets = {'dna', 'rna'};
pattern2.alphawidths = [6 ; 4];
freq2 = cell2mat(arrayfun(@(i) normvec(0.01 + 0.98*rand(5, 1)), (1:sum(pattern2.alphawidths))', 'uniform', false));

posw2 = pwfctdkl(freq2, pattern2);
% invert to find out what background distributions were used
tmp = exp(posw2 - log(freq2./(1 - freq2)));
bkgfreq2 = 1./(1+tmp);
bkgfreq2 = bkgfreq2(:);

bkgfreq2_dna = bkgfreq2(1:4*pattern2.alphawidths(1));
bkgfreq2_rna = bkgfreq2((4*pattern2.alphawidths(1)+1):4*sum(pattern2.alphawidths));

utexpect(all(abs(bkgfreq2_dna - repmat(bkgfreq2_dna(1:4), pattern2.alphawidths(1), 1)) < eps) && ...
    all(abs(bkgfreq2_rna - repmat(bkgfreq2_rna(1:4), pattern2.alphawidths(2), 1)) < eps), ...
    'pwfctdkl multialphabet; dna, rna site consistency');

utexpect(all(bkgfreq2 > 0.1 & bkgfreq2 < 0.5), 'pwfctdkl dna, rna reasonable bkgfreq');

% check treatment of degenerate frequencies -- won't check exactly
% what happens to them, just that we don't get any infinities or nans
freq1(randperm(length(freq1), 23)) = 0;
freq1(randperm(length(freq1), 31)) = 1;
posw1a = pwfctdkl(freq1, pattern1);
freq2(randperm(length(freq2), 6)) = 0;
freq2(randperm(length(freq2), 4)) = 1;
posw2a = pwfctdkl(freq2, pattern2);

utexpect(~any(isnan(posw1a)) && ~any(isnan(posw2a)) && ...
    ~any(isinf(posw1a)) && ~any(isinf(posw2a)), ...
    'pwfctdkl handling degenerate frequencies');

end