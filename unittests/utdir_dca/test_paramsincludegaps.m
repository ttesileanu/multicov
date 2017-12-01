function test_paramsincludegaps
% TEST_PARAMSINCLUDEGAPS Test paramsincludegaps.m.

timer = tic;
pattern.alphabets = {'multi9', 'protein', 'gapbinary'};
pattern.alphawidths = [9 ; 5 ; 4];
binmap = getbinmap(pattern);
m = randn(binmap{end}(end));
mbig = paramsincludegaps(m, pattern);

params.type = 'maxent';
params.alphabets = pattern.alphabets;
params.alphawidths = pattern.alphawidths;
params.refseq(1).seqdb = 'generic';
params.refseq(1).seqid = '';
params.refseq(1).map = 1:9;
params.refseq(2).seqdb = 'generic';
params.refseq(2).seqid = '';
params.refseq(2).map = 1:5;
params.refseq(3).seqdb = 'generic';
params.refseq(3).seqid = '';
params.refseq(3).map = 1:4;

params.couplings = m;
[paramsbig, changed1] = paramsincludegaps(params);

gappattern = alnincludegaps(pattern);

utexpect(isequal(paramsbig.alphabets(:), gappattern.alphabets(:)) && ...
    isequal(paramsbig.alphawidths(:), gappattern.alphawidths(:)) && ...
    all(abs(mbig(:) - paramsbig.couplings(:)) < eps), ...
    'paramsincludegaps same with matrix or params structure', timer);

[~, changed2] = paramsincludegaps(paramsbig);
utexpect(changed1 && ~changed2, 'paramsincludegaps ''changed'' output argument');

gapbinmap = getbinmap(gappattern);
mbigexp = zeros(gapbinmap{end}(end));
n = length(binmap);
for i = 1:n
    idxs1 = binmap{i};
    idxs1big = gapbinmap{i};
    idxs1match = idxs1big(1+(length(idxs1)~=length(idxs1big)):end);
    for j = 1:n
        idxs2 = binmap{j};
        idxs2big = gapbinmap{j};
        idxs2match = idxs2big(1+(length(idxs2)~=length(idxs2big)):end);
        
        mbigexp(idxs1match, idxs2match) = m(idxs1, idxs2);
    end
end
utexpect(all(size(mbig) == gapbinmap{end}(end)) && ...
    all(abs(mbig(:) - mbigexp(:)) < eps), 'paramsincludegaps with random matrix');

end