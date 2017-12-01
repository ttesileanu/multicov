function test_changegauge
% TEST_CHANGEGAUGE Test changegauge.m.

tol = 1e-9;

pattern.alphabets = {'gaprna', 'gapprotein', 'gapdna'};
pattern.alphawidths = [15 ; 7 ; 8];
binmap = getbinmap(pattern);
len = binmap{end}(end);
m = randn(len, len);
m = zeroblockdiag(m + m', pattern);
gauge1 = randn(len, 1);
gauge1cell = cellfun(@(idxs) gauge1(idxs), binmap, 'uniform', false);

params = pattern;
params.type = 'maxent';
params.couplings = m;

params_gauged = changegauge(params, gauge1);
params_gauged_a = changegauge(params, gauge1cell);
utexpect(checkgauge(params_gauged.couplings, params_gauged, gauge1, tol) && ...
    all(abs(params_gauged.couplings(:) - params_gauged_a.couplings(:)) < tol), ...
    'changegauge per-site gauge');

gauge2rna = randn(5, 1);
gauge2prot = randn(21, 1);
gauge2dna = randn(5, 1);
gauge2 = [repmat(gauge2rna, pattern.alphawidths(1), 1) ; ...
          repmat(gauge2prot, pattern.alphawidths(2), 1) ; ...
          repmat(gauge2dna, pattern.alphawidths(3), 1) ...
         ];
m = randn(len, len);
m = zeroblockdiag(m + m', pattern);
params.couplings = m;
params_gauged = changegauge(params, {gauge2rna, gauge2prot, gauge2dna});
utexpect(checkgauge(params_gauged.couplings, params_gauged, gauge2, tol), ...
    'changegauge per-alphabet gauge');

gauge_zs = ones(size(gauge2));
m = randn(len, len);
m = zeroblockdiag(m + m', pattern);
params.couplings = m;
params_gauged = changegauge(params, 'zerosum');
utexpect(checkgauge(params_gauged.couplings, params_gauged, gauge_zs, tol), ...
    'changegauge zero-sum gauge');

gauge_gap = cell2mat(cellfun(@(idxs) [1 ; zeros(length(idxs)-1, 1)], binmap(:), 'uniform', false));
m = randn(len, len);
m = zeroblockdiag(m + m', pattern);
params.couplings = m;
params_gauged = changegauge(params, 'gap');
utexpect(checkgauge(params_gauged.couplings, params_gauged, gauge_gap, tol), ...
    'changegauge gap gauge');

% try single alphabets now
pattern2.alphabets = {'gapmulti5'};
pattern2.alphawidths = 32;
nletts = length(alphagetletters(pattern2.alphabets{1}, 'nogap'));
len2 = pattern2.alphawidths*nletts;

params2 = pattern2;
params2.type = 'maxent';
m2 = randn(len2, len2);
m2 = zeroblockdiag(m2 + m2', pattern2);
params2.couplings = m2;
gauge_1alpha = randn(nletts, 1);
params2_gauged = changegauge(params2, gauge_1alpha);
params2_gauged_a = changegauge(params2, {gauge_1alpha});
utexpect(checkgauge(params2_gauged.couplings, params2_gauged, repmat(gauge_1alpha, params2.alphawidths, 1), tol) && ...
    all(abs(params2_gauged.couplings(:) - params2_gauged_a.couplings(:)) < tol), ...
    'changegauge single alphabet');

% XXX should check that the change is actually a symmetry!

end

function res = checkgauge(m, structure, gaugevec, tol)
% CHECKGAUGE Check that the matrix is in the specified gauge.

res = true;
binmap = getbinmap(structure);
fields = diag(m);
for i = 1:length(binmap)
    idxsi = binmap{i};
    for j = setdiff(1:length(binmap), i)
        idxsj = binmap{j};
        shouldvanish = m(idxsi, idxsj)*flatten(gaugevec(idxsj));
        if any(abs(shouldvanish) > tol)
            res = false;
            break;
        end
    end
    if ~res
        break;
    end
    
    shouldvanish = dot(fields(idxsi), flatten(gaugevec(idxsi)));
    if any(abs(shouldvanish) > tol)
        res = false;
        break;
    end
end

end