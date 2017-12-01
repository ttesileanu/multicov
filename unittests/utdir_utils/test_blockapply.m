function test_blockapply
% TEST_BLOCKAPPLY Test blockapply.m.

% 1. matrices
% first uniform mode
timer = tic;
pattern1.alphabets = {'protein'};
pattern1.alphawidths = 16;
m1 = randn(pattern1.alphawidths*20);
m1 = m1 + m1';
m1red = blockapply(m1, @(m) norm(m, 'fro'), pattern1);
m1red_ns = blockapply(m1, @(m) norm(m, 'fro'), pattern1, 'symmetric', false);
m1red_exp = zeros(pattern1.alphawidths);
for i = 1:pattern1.alphawidths
    idxsi = 20*(i-1) + (1:20);
    for j = 1:pattern1.alphawidths
        idxsj = 20*(j-1) + (1:20);
        m1red_exp(i, j) = norm(m1(idxsi, idxsj), 'fro');
    end
end
% XXX I guess the Frobenius norm of the transpose differs by as much as
% 256*eps from the Frobenius norm of the matrix? why?
utexpect(all(abs(m1red_ns(:) - m1red(:)) < 256*eps) && ...
    all(abs(m1red(:) - m1red_exp(:)) < 256*eps), 'blockapply matrix, single alphabet, symmetric, uniform', timer);

timer = tic;
pattern2.alphabets = {'dna', 'protein', 'rna'};
pattern2.alphawidths = [6 ; 3 ; 8];
m2 = randn(dot(pattern2.alphawidths, [4 ; 20 ; 4]));
m2red = blockapply(m2, @norm, pattern2, 'symmetric', false);
m2red_exp = zeros(size(m2red));
binmap = getbinmap(pattern2);
for i = 1:sum(pattern2.alphawidths)
    for j = 1:sum(pattern2.alphawidths)
        m2red_exp(i, j) = norm(m2(binmap{i}, binmap{j}));
    end
end
utexpect(all(abs(m2red(:) - m2red_exp(:)) < eps), ...
    'blockapply matrix, multi alphabet, not symmetric, uniform', timer);

% now non-uniform -- autodetect, and force
timer = tic;
m1p1 = blockapply(m1, @(m) m + 1, pattern1);
testres = true;
for i = 1:pattern1.alphawidths
    idxsi = 20*(i-1) + (1:20);
    for j = 1:pattern1.alphawidths
        idxsj = 20*(j-1) + (1:20);
        subblock = m1(idxsi, idxsj) + 1;
        if any(abs(m1p1{i, j}(:) - subblock(:)) >= eps)
            testres = false;
            break;
        end
    end
    if ~testres
        break;
    end
end
utexpect(testres, 'blockapply matrix, non-uniform, autodetect', timer);

timer = tic;
m2svals = blockapply(m2, @svd, pattern2, 'symmetric', false, 'uniform', false);
testres = true;
for i = 1:sum(pattern2.alphawidths)
    for j = 1:sum(pattern2.alphawidths)
        subblock = m2(binmap{i}, binmap{j});
        svals = svd(subblock);
        if any(abs(m2svals{i, j} - svals) >= eps)
            testres = false;
            break;
        end
    end
    if ~testres
        break;
    end
end
utexpect(testres, 'blockapply matrix, non-uniform, forced', timer);

% 2. vectors
% uniform
v1 = randn(pattern1.alphawidths*20, 1);
v1red = blockapply(v1, @sum, pattern1);
v1red_exp = zeros(pattern1.alphawidths, 1);
for i = 1:pattern1.alphawidths
    idxsi = 20*(i-1) + (1:20);
    v1red_exp(i) = sum(v1(idxsi));
end
utexpect(all(abs(v1red - v1red_exp) < eps), 'blockapply vector, single alphabet, uniform', timer);

% non-uniform
v2 = randn(dot(pattern2.alphawidths, [4 ; 20 ; 4]), 1);
v2a = blockapply(v2, @flipud, pattern2);
v2a_force = blockapply(v2, @flipud, pattern2, 'uniform', false);
testres = true;
for i = 1:sum(pattern2.alphawidths)
    if ~isequal(v2a{i}, v2a_force{i}) || ~isequal(v2a{i}, flipud(v2(binmap{i})))
        testres = false;
        break;
    end
end
utexpect(testres, 'blockapply vector, multi alphabet, non-uniform', timer);

% 3. indices in applied function
v = randn(binmap{end}(end), 1);
vred_widxs = blockapply(v, @(m, i) i*norm(m), pattern2, 'indices', true);
vexp = zeros(sum(pattern2.alphawidths), 1);
for i = 1:length(vexp)
    vexp(i) = norm(v(binmap{i}))*i;
end
utexpect(all(abs(vexp(:) - vred_widxs(:)) < eps), 'blockapply vector with ''indices'', true');

timer = tic;
m = randn(binmap{end}(end), binmap{end}(end));
m = m + m';
mred_widxs = blockapply(m, @(m, i, j) i*m + j*m.^2, pattern2, 'indices', true, 'symmetric', false);
mexp = cell(sum(pattern2.alphawidths));
testres = true;
for i = 1:size(mexp, 1)
    iidxs = binmap{i};
    for j = 1:size(mexp, 2)
        jidxs = binmap{j};
        block = m(iidxs, jidxs);
        mexp{i, j} = i*block + j*block.^2;
        if any(abs(mexp{i, j}(:) - mred_widxs{i, j}(:)) > eps)
            testres = false;
            break;
        end
    end
    if ~testres
        break;
    end
end
utexpect(testres, 'blockapply matrix with ''indices'' true, non-uniform, non-symmetric', timer);

end