function test_getmaxentenergies
% TEST_GETMAXENTENERGIES Test getmaxentenergies.m.

timer = tic;
tol = 1e-12;

% make up an alignment
N = 64;
pattern.alphabets = {'dna', 'rna', 'protein'};
pattern.alphawidths = [7 ; 11 ; 4];
alignment = alngenrandom(N, pattern);

% make up some couplings
gappattern = alphaincludegaps(pattern);
params = gappattern;
params.type = 'maxent';
params.refseq = struct('seqdb', {'generic', 'generic', 'generic'}, ...
    'seqid', {'', '', ''}, 'map', {1:7, 1:11, 1:4});
binmap = getbinmap(gappattern);
len = binmap{end}(end);
n = length(binmap);

m = randn(len, len);
m = m + m';
params.couplings = zeroblockdiag(m, params) + diag(diag(m));

% get the energies
energies = getmaxentenergies(alignment, params);
binalign = alntobin(alignment);
energies_a = getmaxentenergies(binalign, params);

energies_b = getmaxentenergies(alignment.data, params);

% calculate our own estimates for the energies
energies_exp = zeros(size(energies));
alphabetletters = arrayfun(@(i) repmat({alphagetletters(gappattern.alphabets{i}, 'nogap')}, gappattern.alphawidths(i), 1), ...
    1:length(gappattern.alphabets), 'uniform', false);
alphabetletters = vertcat(alphabetletters{:});
for k = 1:length(energies_exp)
    energy = 0;
    seq = alignment.data(k, :);
    for i = 1:n
        lettersi = alphabetletters{i};
        numi = find(lettersi == seq(i), 1);
        idxi = binmap{i}(numi);
        energy = energy + params.couplings(idxi, idxi)/2;
        for j = (i+1):n
            lettersj = alphabetletters{j};
            numj = find(lettersj == seq(j), 1);
            idxj = binmap{j}(numj);
            energy = energy + params.couplings(idxi, idxj);
        end
    end
    energies_exp(k) = -energy;
end
utexpect(all(abs(energies(:) - energies_exp(:)) < tol), 'getmaxentenergies character alignment', timer);

utexpect(all(abs(energies(:) - energies_a(:)) < tol), 'getmaxentenergies binary alignment');
utexpect(all(abs(energies(:) - energies_b(:)) < tol), 'getmaxentenergies raw sequence data');

% check gauge invariance
timer = tic;
params_zs = changegauge(params, 'zerosum');
params_gap = changegauge(params, 'gap');
energies_zs = getmaxentenergies(alignment, params_zs);
energies_gap = getmaxentenergies(alignment, params_gap);
utexpect(std(energies(:) - energies_zs(:)) < tol && ...
    std(energies(:) - energies_gap(:)) < tol, ...
    'getmaxentenergies gauge invariance', timer);

global multicov_ut_cpppolicy;

if getmaxentenergies('hascpp')
    energies_nocpp = getmaxentenergies(alignment, params, 'nocpp', true);
    energies_b_nocpp = getmaxentenergies(alignment.data, params, 'nocpp', true);
    
    utexpect(all(abs(energies(:) - energies_nocpp(:)) < tol) && ...
        all(abs(energies_b(:) - energies_b_nocpp(:)) < tol), ...
        'getmaxentenergies Matlab code (no C++)');
    
    energies_zs_nocpp = getmaxentenergies(alignment, params_zs, 'nocpp', true);
    energies_gap_nocpp = getmaxentenergies(alignment, params_gap, 'nocpp', true);
    utexpect(std(energies_nocpp(:) - energies_zs_nocpp(:)) < tol && ...
        std(energies_nocpp(:) - energies_gap_nocpp(:)) < tol, ...
        'getmaxentenergies Matlab code gauge invariance', timer);
else
    if strcmp(multicov_ut_cpppolicy, 'warn')
        warning([mfilename ':nocpp'], 'Cannot find C++ extension for getmaxentenergies. Please compile the extension for optimal performance.');
    elseif strcmp(multicov_ut_cpppolicy, 'fail')
        utexpect(false, 'getmaxentenergies C++ extension missing');
    end
end

end