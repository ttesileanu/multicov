function test_setdbbkg
% TEST_SETDBBKG Test setdbbkg.m.

oldfreqs = getdbbkg('multi5');
oldfreqsp = getdbbkg('protein');

setdbbkg('multi5', [0.1 0.2 0.2 0.1 0.4]);
utexpect(isequal(flatten(getdbbkg('multi5')), [0.1 ; 0.2 ; 0.2 ; 0.1 ; 0.4]) && ...
    isequal(getdbbkg('protein'), oldfreqsp), 'setdbbkg new frequencies');

setdbbkg('multi5', []);
utexpect(isequal(getdbbkg('multi5'), oldfreqs) && isequal(getdbbkg('protein'), oldfreqsp), ...
    'setdbbkg reset to default');

end