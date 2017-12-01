function test_paramsexcludegaps
% TEST_PARAMSEXCLUDEGAPS Test paramsexcludegaps.m.

timer = tic;
pattern.alphabets = {'gapmulti9', 'gapprotein', 'binary'};
pattern.alphawidths = [9 ; 5 ; 4];
binmap = getbinmap(pattern);
m = randn(binmap{end}(end));
msmall = paramsexcludegaps(m, pattern);

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
[paramssmall, changed1] = paramsexcludegaps(params);

ungappattern = alnexcludegaps(pattern);

utexpect(isequal(paramssmall.alphabets(:), ungappattern.alphabets(:)) && ...
    isequal(paramssmall.alphawidths(:), ungappattern.alphawidths(:)) && ...
    all(abs(msmall(:) - paramssmall.couplings(:)) < eps), ...
    'paramsexcludegaps same with matrix or params structure', timer);

[~, changed2] = paramsexcludegaps(paramssmall);
utexpect(changed1 && ~changed2, 'paramsexcludegaps ''changed'' output argument');

ungapbinmap = getbinmap(ungappattern);
msmallexp = zeros(ungapbinmap{end}(end));
n = length(binmap);
for i = 1:n
    idxs1 = binmap{i};
    idxs1small = ungapbinmap{i};
    idxs1match = idxs1(1+(length(idxs1)~=length(idxs1small)):end);
    for j = 1:n
        idxs2 = binmap{j};
        idxs2small = ungapbinmap{j};
        idxs2match = idxs2(1+(length(idxs2)~=length(idxs2small)):end);
        
        msmallexp(idxs1small, idxs2small) = m(idxs1match, idxs2match);
    end
end
utexpect(all(size(msmall) == ungapbinmap{end}(end)) && ...
    all(abs(msmall(:) - msmallexp(:)) < eps), 'paramsexcludegaps with random matrix');

end