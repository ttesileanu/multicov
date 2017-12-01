function test_reducecouplings
% TEST_REDUCECOUPLINGS Test reducecouplings.m.

% general-purpose tolerance
tol = 1e-12;

timer = tic;
params.type = 'maxent';
params.alphabets = {'gapdna', 'gapmulti9'};
params.alphawidths = [8 ; 8];
params.refseq = struct('seqdb', {'generic', 'generic'}, 'seqid', {'', ''}, ...
    'map', {1:8, 1:8});

binmap = getbinmap(params);
len = binmap{end}(end);
m = randn(len, len);
m = m + m';
params.couplings = m;
normalize_freq = @(v) v/sum(v);
freq1 = cell2mat(cellfun(@(idxs) normalize_freq(0.02 + 0.96*rand(length(idxs), 1)), ...
    binmap(:), 'uniform', false));

dmodel1 = calculatedirectmodel(params, freq1);
redJ1 = reducecouplings(params, 'di', 'freq1', freq1);

redJ1_exp = zeros(size(redJ1));
n = length(binmap);
for i = 1:n
    idxsi = binmap{i};
    freqsi = freq1(idxsi);
    for j = (i+1):n
        idxsj = binmap{j};
        freqsj = freq1(idxsj);
        block = dmodel1.freq2(idxsi, idxsj);
        prodblock = freqsi(:)*freqsj(:)';
        redJ1_exp(i, j) = trace(block'*log(block./prodblock));
    end
end
redJ1_exp = redJ1_exp + redJ1_exp';
utexpect(all(abs(redJ1(:) - redJ1_exp(:)) < tol) && ...
    all(abs(diag(redJ1)) < eps), 'reducecouplings DI', timer);

redJ1_exp_2 = zeroblockdiag(getmi(dmodel1), 1);
utexpect(all(abs(redJ1_exp_2(:) - redJ1(:)) < tol), ...
    'reducecouplings DI vs. calculatedirectmodel+getmi', timer);

timer = tic;
params_nogap = excludegaps(changegauge(params, 'gap'));
redJ2 = reducecouplings(params, 'norm');
redJ2a = reducecouplings(params_nogap, 'norm');
params_zs = changegauge(params, 'zerosum');
redJ2_exp = zeros(size(redJ2));
for i = 1:n
    idxsi = binmap{i};
    for j = (i+1):n
        idxsj = binmap{j};
        block = params_zs.couplings(idxsi, idxsj);
        redJ2_exp(i, j) = norm(block, 'fro');
    end
end
redJ2_exp = redJ2_exp + redJ2_exp';

utexpect(all(abs(redJ2(:) - redJ2_exp(:)) < tol) && ...
    all(abs(redJ2(:) - redJ2a(:)) < tol) && ...
    all(abs(diag(redJ2)) < eps), 'reducecouplings Frobenius, zero-sum gauge, with/without gaps', timer);

timer = tic;
redJ3 = reducecouplings(params, 'norm', 'gauge', 'original', 'norm', 'spectral');
redJ3a = reducecouplings(params, 'norm', 'gauge', 'original', 'norm', @norm);
redJ3_exp = zeros(size(redJ3));
for i = 1:n
    idxsi = binmap{i};
    for j = (i+1):n
        idxsj = binmap{j};
        block = params.couplings(idxsi, idxsj);
        redJ3_exp(i, j) = norm(block);
    end
end
redJ3_exp = redJ3_exp + redJ3_exp';

utexpect(all(abs(redJ3(:) - redJ3_exp(:)) < tol) && ...
    all(abs(redJ3(:) - redJ3a(:)) < eps) && ...
    all(abs(diag(redJ3)) < eps), 'reducecouplings spectral, original gauge', timer);

end