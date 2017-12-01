function test_findfiles
% TEST_FINDFILES Test findfiles.m.

global multicov_ut_dir;

cd(multicov_ut_dir);

files = findfiles;
utexpect(any(ismember(files, fullfile('utdir_utils', 'test_findfiles.m'))), 'findfiles');

end