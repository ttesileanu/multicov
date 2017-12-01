function test_seqtruncate
% TEST_SEQTRUNCATE Test seqtruncate.m.

alignment1 = alnmake(['AA-U' ; 'ACGU' ; '--GG' ; 'A-GG' ; 'GGC-'], 'dna');
newalign1 = seqtruncate(alignment1);
newalign1_exp = alnmake(['AAU' ; 'ACU' ; '--G' ; 'A-G' ; 'GG-'], 'dna');
%newalign1_exp.refseq.seq = 'AAU';
utexpect(isequal(newalign1, newalign1_exp), 'seqtruncate to first sequence');

newalign2 = seqtruncate(alignment1, 3, 'truncate', false);
refseq_exp.seqdb = 'generic';
refseq_exp.seqid = '';
%refseq_exp.seq = 'GG';
refseq_exp.map = [0 ; 0 ; 1 ; 2];
utexpect(isequal(newalign2.refseq, refseq_exp) && ...
    isequal(newalign2.data, alignment1.data), ...
    'seqtruncate to arbitrary sequence, w/o truncation');

newalign3 = seqtruncate(alignment1, 2, 'seq', 'GUAACCGUU');
refseq_exp.seqdb = 'generic';
refseq_exp.seqid = '';
%refseq_exp.seq = 'GUAACCGUU';
refseq_exp.map = [4 ; 6 ; 7 ; 8];
utexpect(isequal(newalign3.refseq, refseq_exp), 'seqtruncate to given sequence');

alignment2.alphabets = {'protein', 'rna'};
alignment2.alphawidths = [6 ; 5];
alignment2.data = [...
    '-AWGGHGG-CA' ; ...
    'D-GG-AC-ACC' ; ...
    'WWGYPDGAACC' ; ...
    'W--IIK-UUAU' ; ...
    '--FDGHUCGAU' ; ...
];
alignment2.refseq(1).seqdb = 'generic';
alignment2.refseq(2).seqdb = 'generic';
alignment2.refseq(1).seqid = '';
alignment2.refseq(2).seqid = '';
alignment2.refseq(1).map = (1:6)';
alignment2.refseq(2).map = (1:5)';
alignment2.seqw = ones(5, 1);
alignment2.annotations = {'' ; '' ; '' ; '' ; ''};
alignment2.type = 'character';
expseq1 = 'AWCWGYPPDY';
expseq2 = 'AUGACACGCUU';
[newalign4, alignfrac, aligndata] = seqtruncate(alignment2, 3, 'seq', {expseq1, expseq2});
refseq_exp(1).seqdb = 'generic';
refseq_exp(2).seqdb = 'generic';
refseq_exp(1).seqid = '';
refseq_exp(2).seqid = '';
refseq_exp(1).map = [2 ; 4 ; 5 ; 6 ; 8 ; 9];
refseq_exp(2).map = [3 ; 4 ; 6 ; 7 ; 9];

[~, align_prot] = nwalign(alignment2.data(3, 1:6), expseq1, 'glocal', true);
[~, align_rna] = nwalign(alignment2.data(3, 7:end), expseq2, 'glocal', true, 'alphabet', 'nt');

gapalignment2 = alnincludegaps(alignment2);
gapnewalign4 = seqtruncate(gapalignment2, 3, 'seq', {expseq1, expseq2});
gapnewalign4_exp = alnincludegaps(newalign4);
utexpect(isequal(newalign4.refseq, refseq_exp) && ...
    isequal(newalign4.data, alignment2.data) && ...
    strcmp(aligndata{1}.alignment, align_prot) && ...
    strcmp(aligndata{2}.alignment, align_rna) && ...
    abs(alignfrac - 1) < eps && ...
    isequal(gapnewalign4_exp, gapnewalign4), ...
    'seqtruncate with multiple alphabets, and alignfrac, aligndata, with/without gapped alphabets');

newalign4a0 = seqtruncate(alignment2, 1, 'seq', {'AGHAAC', 'GCAACTCT'}, 'verbose', false);
utexpect(isequal(newalign4a0.data, alignment2.data(:, [4 5 6 8 10 11])), ...
    'seqtruncate check data with multiple alphabets');

[~, alignfrac, ~] = seqtruncate(alignment2, 3, 'seq', {'AWCWGYPPCY', 'AUGACACGCUU'}, ...
    'verbose', false);
utexpect(abs(alignfrac - 10/11) < eps, 'seqtruncate with imperfect align');

newalign3a = seqtruncate(alignment1, 2, 'seq', 'GUAACCGUU', ...
    'posnames', {'1', '2', '3', '3a', '4', '5', '6', '8', '9'});
refseq3a_exp.seqdb = 'generic';
refseq3a_exp.seqid = '';
%refseq3a_exp.seq = 'GUAACCGUU';
refseq3a_exp.map = {'3a' ; '5' ; '6' ; '8'};
utexpect(isequal(newalign3a.refseq, refseq3a_exp), 'seqtruncate with posnames');

newalign4a = seqtruncate(alignment2, 3, 'seq', {'AWCWGYPPDY', 'AUGACACGCUU'}, ...
    'posnames', {[], [2 3 4 5 7 8 9 10 11 13 14]});
refseq_exp(1).seqdb = 'generic';
refseq_exp(2).seqdb = 'generic';
refseq_exp(1).seqid = '';
refseq_exp(2).seqid = '';
refseq_exp(1).map = [2 ; 4 ; 5 ; 6 ; 8 ; 9];
refseq_exp(2).map = [4 ; 5 ; 8 ; 9 ; 11];
%refseq_exp(1).seq = 'AWCWGYPPDY';
%refseq_exp(2).seq = 'AUGACACGCUU';
utexpect(isequal(newalign4a.refseq, refseq_exp), 'seqtruncate with posnames, multiple alphabets');

end