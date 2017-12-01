function test_tofieldname
% TEST_TOFIELDNAME Test tofieldname.m.

utexpect(strcmp(tofieldname('foo'), 'foo'), 'tofieldname simple');
utexpect(strcmp(tofieldname('foo-bar'), 'foo_bar'), 'tofieldname dashes');
utexpect(strcmp(tofieldname('Who Is This?'), 'Who_Is_This_'), 'tofieldname spaces, symbols');
utexpect(strcmp(tofieldname('-_foo#'), 'foo_'), 'tofieldname start with letter');

end