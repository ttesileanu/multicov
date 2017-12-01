function test_calculatedirectmodel
% TEST_CALCULATEDIRECTMODEL Test calculatedirectmodel.m.

% general-use tolerance
tol = 1e-12;
% tolerance used for marginals
mtol = 1e-3;

timer = tic;
params.type = 'maxent';
params.alphabets = {'gaprna', 'gapprotein'};
params.alphawidths = [9 ; 7];
params.refseq = struct('seqdb', {'generic', 'generic'}, 'seqid', {'', ''}, ...
    'map', {1:9, 1:7});

binmap = getbinmap(params);
len = binmap{end}(end);
m = randn(len, len);
m = m + m';
params.couplings = m;
normalize_freq = @(v) v/sum(v);
freq1 = cell2mat(cellfun(@(idxs) normalize_freq(0.02 + 0.96*rand(length(idxs), 1)), ...
    binmap(:), 'uniform', false));

dmodel1 = calculatedirectmodel(params, freq1);
% test that the block-diagonal is empty and that the marginals are correct
utexpect(checkmodel(dmodel1, tol, mtol) && all(isfield(dmodel1, {'freq1', 'freq2', 'alphabets', 'alphawidths', 'refseq'})) && ...
    isequal(dmodel1.alphabets(:), params.alphabets(:)) && ...
    isequal(dmodel1.alphawidths(:), params.alphawidths(:)) && ...
    isequal(dmodel1.refseq(:), params.refseq(:)) && ...
    isequal(dmodel1.freq1(:), freq1(:)), ...
    'calculatedirectmodel basic test', timer);

% test that the block-diagonal elements are inconsequential
timer = tic;
params.couplings = zeroblockdiag(params.couplings, params);
dmodel1a = calculatedirectmodel(params, freq1);
comparemodels = @(m1, m2, tol1, tol2) ...
    all(abs(m1.freq1(:) - m2.freq1(:)) < tol1) && ...
    all(abs(m1.freq2(:) - m2.freq2(:)) < tol2);
utexpect(comparemodels(dmodel1, dmodel1a, tol, tol), ...
    'calculatedirectmodel independence on maxent model fields', timer);

% test that this works when some of the inputs don't include gaps
timer = tic;
params_gapgauge = changegauge(params, 'gap');
params_nogap = paramsexcludegaps(params_gapgauge);
freq1_nogap = statsexcludegaps(freq1, params);

dmodel2 = calculatedirectmodel(params_nogap, freq1);
dmodel3 = calculatedirectmodel(params_nogap, freq1_nogap);
dmodel4 = calculatedirectmodel(params, freq1_nogap);
utexpect(comparemodels(dmodel2, dmodel3, tol, mtol) && comparemodels(dmodel2, dmodel4, tol, mtol) && ...
    comparemodels(dmodel1, dmodel4, tol, mtol), ...
    'calculatedirectmodel gapless inputs and gauge invariance', timer);

% check tol input
timer = tic;
dmodel5 = calculatedirectmodel(params, freq1, 'tol', 1e-6);
utexpect(checkmodel(dmodel5, tol, 1e-5), ...
    'calculatedirectmodel different ''tol''', timer);

% check nfailed, tol, and maxiter inputs
% to avoid randomly-occurring problems, fix the random seed
timer = tic;
rng('default');
m = randn(len, len);
m = m + m';
params.couplings = m;
freq1 = cell2mat(cellfun(@(idxs) normalize_freq(0.02 + 0.96*rand(length(idxs), 1)), ...
    binmap(:), 'uniform', false));

[~, nfailed] = calculatedirectmodel(params, freq1, 'tol', 1e-6, 'maxiter', 10, 'verbose', false);
utexpect(nfailed > 0, 'calculatedirectmodel ''nfailed'' output, and ''maxiter'' option', timer);

% XXX need to check nfailed output argument, tol and maxiter inputs


end

function testres = checkmodel(dmodel, tol, mtol)
% CHECKMODEL Check the maxent model.

binmap = getbinmap(dmodel);
freq1 = dmodel.freq1;
n = length(binmap);
testres = true;
for i = 1:n
    idxsi = binmap{i};
    block = dmodel.freq2(idxsi, idxsi);
    if any(block(:) > tol)
        testres = false;
        break;
    end
    for j = setdiff(1:n, i)
        idxsj = binmap{j};
        block = dmodel.freq2(idxsi, idxsj);
        Pi = sum(block, 2);
        Pj = sum(block, 1);
        if any(abs(Pi(:) - flatten(freq1(idxsi))) > mtol) || any(abs(Pj(:) - flatten(freq1(idxsj))) > mtol)
            testres = false;
            break;
        end
    end
    if ~testres
        break;
    end
end

end