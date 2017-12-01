function test_guesstype
% TEST_GUESSTYPE Test guesstype.m.

[t1, v1] = guesstype('');
utexpect(strcmp(t1, 'x') && isempty(v1), 'guesstype x (empty)');

[t2, v2] = guesstype('tRue');
[t2a, v2a] = guesstype('false');
[t2b, v2b] = guesstype('False');
utexpect(strcmp(t2, 'l') && strcmp(t2a, 'l') && strcmp(t2b, 'l') && ...
    islogical(v2) && islogical(v2a) && islogical(v2b) && ...
    isscalar(v2) && isscalar(v2a) && isscalar(v2b) && ...
    v2 && ~v2a && ~v2b, 'guesstype l (logical)');

[t3, v3] = guesstype('5');
[t3a, v3a] = guesstype('3.14');
utexpect(strcmp(t3, 'n') && strcmp(t3a, 'n') && abs(v3 - 5) < eps && abs(v3a - 3.14) < eps, ...
    'guesstype n (numeric scalar)');

[t3b, v3b] = guesstype('nAn');
utexpect(strcmp(t3b, 'n') && isnan(v3b), 'guesstype nan (scalar)');

[t4, v4] = guesstype('3 5 7 10');
utexpect(strcmp(t4, 'v') && isequal(v4, [3 5 7 10]), ...
    'guesstype v (numeric vector');

[t4a, v4a] = guesstype('3 5 nan');
utexpect(strcmp(t4a, 'v') && isequal(v4a(1:2), [3 5]) && isnan(v4a(3)), ...
    'guesstype nan in vector');

[t5, v5] = guesstype('foo');
utexpect(strcmp(t5, 's') && strcmp(v5, 'foo'), ...
    'guesstype s (string)');

end