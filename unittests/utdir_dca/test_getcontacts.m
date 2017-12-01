function test_getcontacts
% TEST_GETCONTACTS Test getcontacts.m.

% test that we're finding contacts in a mock matrix
mockJred = [...
%   1   3   4   6   a   b   c
    2  1.5 0.2 1.2 0.1  0  0.4 ; ...    %1
   1.5  0  0.9 0.1 0.1 0.1 0.2 ; ...    %3
   0.2 0.9  0   0   0  0.2  0  ; ...    %4
   1.2 0.1  0   3  1.4 0.1 0.1 ; ...    %6
   0.1 0.1  0  1.4  0  0.8 0.3 ; ...    %a
    0  0.1 0.2 0.1 0.8  0  1.0 ; ...    %b
   0.4 0.2  0  0.1 0.3 1.0  0    ...    %c
];
structure1.alphabets = {'protein' ; 'dna'};
structure1.alphawidths = [4 ; 3];
structure1.refseq(1).seqdb = 'generic';
structure1.refseq(1).seqid = 'fake';
structure1.refseq(1).map = [1 ; 3 ; 4 ; 6];
structure1.refseq(2).seqdb = 'generic';
structure1.refseq(2).seqid = 'fake';
structure1.refseq(2).map = {'a' ; 'b' ; 'c'};

contacts1 = getcontacts(mockJred, structure1, 'seqcutoff', 0);
utexpect(all(isfield(contacts1, {'mapped', 'matrix', 'rawpairs', 'redJ', 'refseq'})) && ...
    isequal(contacts1.redJ, mockJred) && isequal(contacts1.refseq, structure1.refseq), ...
    'getcontacts output has the right fields; redJ, and refseq are correct');

mat_exp = (mockJred >= 0.4);
mat_exp(logical(eye(size(mat_exp)))) = false;
utexpect(isequal(contacts1.matrix, mat_exp), ...
    'getcontacts ''matrix'' is correct (no seqcutoff)');

utexpect(isequal(contacts1.rawpairs, ...
   [1 4 1 6 2 5 1 ; ...
    2 5 4 7 3 6 7]), 'getcontacts check ''rawpairs'', no seqcutoff');

utexpect(isequal(size(contacts1.mapped), [2 2]) && ...
    isequal(contacts1.mapped{1, 1}, {1 1 3 ; 3 6 4}) && ...
    isequal(contacts1.mapped{1, 2}, {6 1 ; 'a' 'c'}) && ...
    isequal(contacts1.mapped{2, 1}, {'a' 'c' ; 6 1}) && ...
    isequal(contacts1.mapped{2, 2}, {'b' 'a' ; 'c' 'b'}), ...
    'getcontacts check ''mapped'', no seqcutoff');

contacts2 = getcontacts(mockJred, structure1, 'n', 4);
utexpect(isequal(contacts2.matrix, [...
    false false false true  false false false ; ...
    false false false false false false false ; ...
    false false false false false false false ; ...
    true  false false false true  false false ; ...
    false false false true  false true  false ; ...
    false false false false true  false true  ; ...
    false false false false false true  false]), ...
    'getcontacts std. seqcutoff, scalar ''n''');

structure2 = structure1;
structure2.refseq(2).map = [1 ; 3 ; 4];

contacts3 = getcontacts(mockJred, structure2, 'seqcutoff', [2 ; 1], 'n', [2 2 ; 2 1]);
utexpect(isequal(contacts3.mapped{1, 1}, {1 1 ; 6 4}) && ...
    isequal(contacts3.mapped{1, 2}, {6 1 ; 1 4}) && ...
    isequal(contacts3.mapped{2, 1}, {1 4 ; 6 1}) && ...
    isequal(contacts3.mapped{2, 2}, {1 ; 3}), ...
    'getcontacts per-alphabet seqcutoff, per pair of alphabets ''n''');

end