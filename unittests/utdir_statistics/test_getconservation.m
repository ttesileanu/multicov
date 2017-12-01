function test_getconservation
% TEST_GETCONSERVATION Test getconservation.m.

timer = tic;
alignment1 = alngenrandom(1024, 128, 'protein');
binalign1 = alntobin(alignment1);
stats.alphabets = alignment1.alphabets;
stats.alphawidths = alignment1.alphawidths;
stats.freq1 = getfreq1(alignment1);
cons1_aln = getconservation(alignment1, 'method', 'max');
cons1_bin = getconservation(binalign1, 'method', 'max');
cons1_sts = getconservation(stats, 'method', 'max');
cons1_exp = arrayfun(@(i) max(stats.freq1((1+20*(i-1)):20*i)), 1:alignment1.alphawidths);
utexpect(all(abs(cons1_aln - cons1_bin) < eps) && ...
    all(abs(cons1_aln - cons1_sts) < eps) && ...
    all(abs(cons1_aln(:) - cons1_exp(:)) < eps), ...
    'getconservation single alphabet, ''method'' max', timer);

cons1_wg = getconservation(binalign1, 'method', 'max', 'gaps', true);
exp_freq = @(v) [1-sum(v(:)) ; v(:)];
cons1_wg_exp = arrayfun(@(i) max(exp_freq(stats.freq1((1+20*(i-1)):20*i))), 1:alignment1.alphawidths);
utexpect(all(abs(cons1_wg(:) - cons1_wg_exp(:)) < eps), ...
    'getconservation ''method'' max with gaps');

% the harder part: testing the KL divergence code
pattern2.alphabets = {'dna', 'protein', 'rna'};
pattern2.alphawidths = [32 ; 16 ; 24];
alignment2 = alngenrandom(1600, pattern2);
freq2 = getfreq1(alignment2);
cons2_nogap = getconservation(alignment2, 'method', 'kl', 'gaps', false);
cons2_nogap_exp = getklcons(freq2, pattern2, []);
utexpect(all(abs(cons2_nogap(:) - cons2_nogap_exp(:)) < eps), ...
    'getconservation ''method'' kl, no gaps');

gaps = getgapstructure(alignment2);

cons2_autogap = getconservation(alignment2, 'method', 'kl', 'gaps', true, 'bkggap', 'auto');
gapsalpha = [...
    repmat(mean(flatten(gaps(:, 1:32))), 32, 1) ; ...
    repmat(mean(flatten(gaps(:, 33:48))), 16, 1) ; ...
    repmat(mean(flatten(gaps(:, 49:end))), 24, 1) ...
];
cons2_autogap_exp = getklcons(freq2, pattern2, gapsalpha);
utexpect(all(abs(cons2_autogap - cons2_autogap_exp) < eps), ...
    'getconservation ''method'' kl, gaps auto');

cons2_autoallgap = getconservation(alignment2, 'method', 'kl', 'gaps', true, 'bkggap', 'autoall');
gapsall = repmat(mean(gaps(:)), sum(pattern2.alphawidths), 1);
cons2_autoallgap_exp = getklcons(freq2, pattern2, gapsall);
utexpect(all(abs(cons2_autoallgap - cons2_autoallgap_exp) < eps), ...
    'getconservation ''method'' kl, gaps autoall');

gapscalar = 0.05;
cons2_gapscalar = getconservation(alignment2, 'method', 'kl', 'gaps', true, 'bkggap', gapscalar);
cons2_gapscalar_exp = getklcons(freq2, pattern2, repmat(gapscalar, sum(pattern2.alphawidths), 1));
utexpect(all(abs(cons2_gapscalar - cons2_gapscalar_exp) < eps), ...
    'getconservation ''method'' kl, gaps scalar');

gapvecalpha = [0.1 0.05 0.1];
cons2_gapvecalpha = getconservation(alignment2, 'method', 'kl', 'gaps', true, 'bkggap', gapvecalpha);
cons2_gapvecalpha_exp = getklcons(freq2, pattern2, [...
    repmat(gapvecalpha(1), pattern2.alphawidths(1), 1) ; ...
    repmat(gapvecalpha(2), pattern2.alphawidths(2), 1) ; ...
    repmat(gapvecalpha(3), pattern2.alphawidths(3), 1) ...
]);
utexpect(all(abs(cons2_gapvecalpha - cons2_gapvecalpha_exp) < eps), ...
    'getconservation ''method'' kl, gaps per alphabet');

gapvecall = 0.1*rand(sum(alignment2.alphawidths), 1);
cons2_gapvecall = getconservation(alignment2, 'method', 'kl', 'gaps', true, 'bkggap', gapvecall);
cons2_gapvecall_exp = getklcons(freq2, pattern2, gapvecall);
utexpect(all(abs(cons2_gapvecall - cons2_gapvecall_exp) < eps), ...
    'getconservation ''method'' kl, gaps per residue');

timer = tic;
gapalignment1 = alnincludegaps(alignment1);
gapalignment2 = alnincludegaps(alignment2);
gapcons1_aln = getconservation(gapalignment1, 'method', 'max');
gapcons1_wg = getconservation(gapalignment1, 'method', 'max', 'gaps', true);
gapcons2_nogap = getconservation(gapalignment2, 'method', 'kl', 'gaps', false);
gapcons2_autogap = getconservation(gapalignment2, 'method', 'kl', 'gaps', true, 'bkggap', 'auto');
gapcons2_autoallgap = getconservation(gapalignment2, 'method', 'kl', 'gaps', true, 'bkggap', 'autoall');
gapcons2_gapscalar = getconservation(gapalignment2, 'method', 'kl', 'gaps', true, 'bkggap', gapscalar);
gapcons2_gapvecalpha = getconservation(gapalignment2, 'method', 'kl', 'gaps', true, 'bkggap', gapvecalpha);
gapcons2_gapvecall = getconservation(gapalignment2, 'method', 'kl', 'gaps', true, 'bkggap', gapvecall);
testres = cellfun(@(a, b) all(abs(a(:) - b(:)) < eps), ...
    {cons1_aln, cons1_wg, cons2_nogap, cons2_autogap, cons2_autoallgap, cons2_gapscalar, cons2_gapvecalpha, cons2_gapvecall}, ...
    {gapcons1_aln, gapcons1_wg, gapcons2_nogap, gapcons2_autogap, gapcons2_autoallgap, gapcons2_gapscalar, gapcons2_gapvecalpha, gapcons2_gapvecall});
utexpect(all(testres), 'getconservation with gapped alignments', timer);

end

function res = getklcons(freq, pattern, gaps)
% GETKLCONS Calculate KL conservation using the SCA background distribution.
%   XXX This doesn't behave the same as getconservation for gapped alignments!

bkgfreq_prot = [.073 .025 .050 .061 .042 .072 .023 .053 .064 .089 .023 .043 .052 .040 .052 .073 .056 .063 .013 .033];

n = sum(pattern.alphawidths);
binmap = getbinmap(pattern);
res = zeros(n, 1);
for i = 1:n
    alpha = pattern.alphabets{find(cumsum(pattern.alphawidths) >= i, 1)};
    nletters = length(alphagetletters(alpha, 'nogap'));
    if strcmp(alpha, 'protein')
        bkgfreq = bkgfreq_prot;
    else
        bkgfreq = ones(nletters, 1) / nletters;
    end
    
    subfreq = freq(binmap{i});
    if isempty(gaps)
        s = sum(subfreq);
        if abs(s) > eps
            subfreq = subfreq / s;
        end
        subfreq = subfreq(:);
        bkgfreq = bkgfreq(:);
    else
        subfreq = [1 - sum(subfreq(:)) ; subfreq(:)];
        bkgfreq = [gaps(i) ; (1 - gaps(i))*bkgfreq(:)];
    end
    mask = (abs(subfreq) > eps);
    res(i) = sum(subfreq(mask).*log(subfreq(mask) ./ bkgfreq(mask)));
end

end