function test_findcommonidxs
% TEST_FINDCOMMONIDXS Test findcommonidxs.m.

refa_1 = [1 ; 2 ; 3 ; 5 ; 7 ; 9 ; 10 ; 11];
refa_2 = [3 ; 4 ; 5 ; 6 ; 7 ; 8];
[idxsa_1, idxsa_2] = findcommonidxs(refa_1, refa_2);
idxsa_1_exp = [3 ; 4 ; 5];
idxsa_2_exp = [1 ; 3 ; 5];
utexpect(isequal(idxsa_1(:), idxsa_1_exp(:)) && isequal(idxsa_2(:), idxsa_2_exp(:)), ...
    'findcommonidxs numeric vectors');

refb_1 = {'A' ; 'B' ; 'Ex' ; 'F' ; 'Z'};
refb_2 = {'B' ; 'C' ; 'Ex' ; 'W' ; 'Y'};
[idxsb_1, idxsb_2] = findcommonidxs(refb_1, refb_2);
idxsb_1_exp = [2 ; 3];
idxsb_2_exp = [1 ; 3];
utexpect(isequal(idxsb_1(:), idxsb_1_exp(:)) && isequal(idxsb_2(:), idxsb_2_exp(:)), ...
    'findcommonidxs cell array of strings');

[idxsc_1, idxsc_2] = findcommonidxs(refa_1, refb_2);
[idxsd_1, idxsd_2] = findcommonidxs(refb_1, refa_2);
utexpect(isempty(idxsc_1) && isempty(idxsc_2) && isempty(idxsd_1) && isempty(idxsd_2), ...
    'findcommonidxs no overlap between numbers and strings');

[idxsd_1, idxsd_2] = findcommonidxs([1 2 3 0 0 4 5], [2 3 4 0 0 6]);
utexpect(isequal(idxsd_1(:), [2 ; 3 ; 6]) && isequal(idxsd_2(:), [1 ; 2 ; 3]), ...
    'findcommonidxs 0 is not common');

% check that the ordering doesn't confuse findcommonidxs
ref1_disord = [1 ; 5 ; 3 ; 7 ; 2 ; 9 ; 11 ; 10];
ref2_disord = [3 ; 4 ; 6 ; 5 ; 7 ; 8];
[idxs1_disord, idxs2_disord] = findcommonidxs(ref1_disord, ref2_disord);
utexpect(isequal(ref1_disord(idxs1_disord), ref2_disord(idxs2_disord)), ...
    'findcommonidxs when ordering is different between the two sequences');

end