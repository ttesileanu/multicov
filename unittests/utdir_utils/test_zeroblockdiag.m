function test_zeroblockdiag
% TEST_ZEROBLOCKDIAG Test zeroblockdiag.m.

n = 32;
pattern1.alphabets = {'dna'};
pattern1.alphawidths = n;
m = randn(4*n, 4*n);
m_zd = zeroblockdiag(m, pattern1);
m_zd_2 = zeroblockdiag(m, 4);
m_exp = m;
for i = 1:n
    idxs = (i-1)*4 + (1:4);
    m_exp(idxs, idxs) = 0;
end
utexpect(all(abs(m_zd(:) - m_exp(:)) < eps) && ...
    all(abs(m_zd_2(:) - m_exp(:)) < eps), 'zeroblockdiag single alphabet or given block size');

pattern2.alphabets = {'protein', 'rna'};
pattern2.alphawidths = [16 ; 13];
binmap = getbinmap(pattern2);
m = randn(binmap{end}(end), binmap{end}(end));
clearval = 0.5;
m_cleard = zeroblockdiag(m, pattern2, 3, 'value', clearval);
m_exp = m;
n2 = length(binmap);
for i = 1:n2
    idxs1 = binmap{i};
    m_exp(idxs1, idxs1) = clearval;
    for j = (i+1):min(n2, i+3)
        idxs2 = binmap{j};
        m_exp(idxs1, idxs2) = clearval;
        m_exp(idxs2, idxs1) = clearval;
    end
end
utexpect(all(abs(m_cleard(:) - m_exp(:)) < eps), ...
    'zeroblockdiag multi alphabet, with several diagonals, non-zero fill value');

end