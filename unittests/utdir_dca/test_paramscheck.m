function test_paramscheck
% TEST_PARAMSCHECK Test paramscheck.m.

timer = tic;
N = 128;
pattern.alphabets = {'protein', 'multi5', 'binary'};
pattern.alphawidths = [7 ; 10 ; 12];
alignment = alngenrandom(N, pattern);
dca = getdca(alignment);
params = getmaxentparams(dca, 'extended', false);
utexpect(paramscheck(params), 'paramscheck with input from getmaxentparams on random alignment', ...
    timer);

mockp.type = 'maxent';
mockp.alphabets = pattern.alphabets;
mockp.alphawidths = pattern.alphawidths;
mockp.refseq(1).seqdb = 'generic';
mockp.refseq(1).seqid = '';
mockp.refseq(1).map = 1:7;
mockp.refseq(2).seqdb = 'multi5';
mockp.refseq(2).seqid = '';
mockp.refseq(2).map = 1:10;
mockp.refseq(3).seqdb = 'binary';
mockp.refseq(3).seqid = '';
mockp.refseq(3).map = 1:12;
sz = 20*7 + 4*10 + 12;
mockp.couplings = randn(sz, sz);
utexpect(paramscheck(mockp), 'paramscheck with mock structure');

end