function test_matautocorr
% TEST_MATAUTOCORR Test matautocorr.m.

N = 128;
n = 32;

% constant matrix should give perfect autocorrelation
m0 = ones(N, n);
ac0 = matautocorr(m0);
utexpect(all(abs(ac0(:) - ones(size(ac0))) < eps) && ...
    isvector(ac0) && length(ac0) == N/2-1, ...
    'matautocorr with constant matrix');

m1 = repmat(randn(1, n), N, 1);
ac1 = matautocorr(m1);
utexpect(all(abs(ac1(:) - ones(size(ac1))) < eps), ...
    'matautocorr with constant vector');

m2 = randn(N, n);
ac2 = matautocorr(m2, [],  'slow', false);
ac2s = matautocorr(m2, [], 'slow', true);
utexpect(all(abs(ac2(:) - ac2s(:)) < 1e-12), ...
    'matautocorr slow vs. fast version, random matrix');

range = (3:5:N/3);
ac2r = matautocorr(m2, range);
utexpect(all(abs(flatten(ac2(range)) - ac2r(:)) < eps), ...
    'matautocorr with non-default range');

N2 = 512;

m3 = randn(N2, n);
m3(2:end, :) = 0.5*m3(2:end, :);
m3 = cumsum(m3, 1);
ac3 = matautocorr(m3, [], 'slow', false);
ac3s = matautocorr(m3, [], 'slow', true);
utexpect(all(abs(ac3(:) - ac3s(:)) < 1e-12), ...
    'matautocorr slow vs. fast version, random walk');

m4 = zeros(N2, n);
m4(1, :) = randn(1, n);
for i = 2:N2
    m4(i, :) = 0.8*m4(i - 1, :) + 0.1*randn(1, n);
end
range = 5:3:(N/4);
ac4 = matautocorr(m4, range, 'slow', false);
ac4s = matautocorr(m4, range, 'slow', true);
utexpect(all(abs(ac4(:) - ac4s(:)) < 1e-12), ...
    'matautocorr slow vs. fast version, autoregressive');

utexpect(abs(ac4(end)) < 0.1, 'matautocorr small at end of autoregressive process');

period = 48;
m5a = bsxfun(@times, flatten(sin((1:N2)*(2*pi/period))), randn(1, n));
m5b = bsxfun(@times, flatten(cos((1:N2)*(2*pi/period))), randn(1, n));
m5 = m5a + m5b;
ac5 = matautocorr(m5);

ac5fft = fft(ac5);
[~, fftmax] = max(abs(ac5fft(1:floor(end/2))));
utexpect(abs(fftmax - length(ac5)/period) <= 2, 'matautocorr periodic data');

end