function test_splitoptions
% TEST_SPLITOPTIONS Test splitoptions.m.

options = {'a', 3, 'b', 5, 'c', 7, 'alpha', 'x', 'beta', 1, 'gamma', 0};
names1 = {'a', 'c', 'alpha'};
names2 = {'b', 'alpha', 'beta'};
[opts1, opts2, others] = splitoptions(options, names1, names2);
utexpect(isequal(opts1(:), {'a' ; 3 ; 'c' ; 7 ; 'alpha' ; 'x'}) && ...
    isequal(opts2(:), {'b' ; 5 ; 'alpha' ; 'x' ; 'beta' ; 1}) && ...
    isequal(others(:), {'gamma' ; 0}), ...
    'splitoptions');

[opts, others] = splitoptions({}, {'aleph'});
utexpect(isempty(opts) && isempty(others), 'splitoptions empty cell array');

end