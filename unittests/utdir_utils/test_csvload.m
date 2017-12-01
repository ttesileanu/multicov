function test_csvload
% TEST_CSVLOAD Test csvload.m.

global multicov_ut_dir;

fname = fullfile(multicov_ut_dir, 'temp', 'csvload_test.csv');
f = fopen(fname, 'wt');
fprintf(f, 'n,Xx\n');
fprintf(f, '1,true\n');
fprintf(f, '# this is a comment\n');
fprintf(f, ' 2, false\n');
fclose(f);
data = csvload(fname);
utexpect(isequal(data.header(:), {'n' ; 'Xx'}) && ...
    isequal(data.columns(:), {[1 ; 2] ; [true ; false]}) && ...
    isequal(data.coltypes, 'nl') && ...
    all(isfield(data.data, {'n', 'xx'})) && ...
    isequal(data.data.n, [1 ; 2]) && isequal(data.data.xx, [true ; false]), ...
    'csvload autodetect header; n, l columns; trim ws; name processing; comment');

data_raw = csvload(fname, 'process', false);
utexpect(isempty(data_raw.header) && isequal(data_raw.columns(:), {{'n' ; '1' ; '2'} ; {'Xx' ; 'true' ; 'false'}}), ...
    'csvload no processing');

f = fopen(fname, 'w');
fprintf(f, 'foo bar Name+2 empTy\n');
fprintf(f, '1 true text1 \n');
fprintf(f, '%% this is a comment\n');
fprintf(f, '2 0 text2 \n');
fprintf(f, ' true text3 \n');
fprintf(f, '# this is also a comment\n');
fclose(f);
data = csvload(fname, 'delimiter', ' ', 'comment', '%#', 'header', true, ...
    'namefct', []);
utexpect(isequal(data.header(:), {'foo' ; 'bar' ; 'Name+2' ; 'empTy'}) && ...
    isequal(flatten(data.columns(2:end)), ...
        {{'true' ; '0' ; 'true'} ; {'text1' ; 'text2' ; 'text3'} ; {'' ; '' ; ''}}) && ...
    isequal(flatten(data.columns{1}(1:2)), [1 ; 2]) && isnan(data.columns{1}(3)) && ...
    isequal(data.coltypes, 'nssx') && ...
    all(isfield(data.data, {'foo', 'bar', 'Name_2', 'empTy'})) && ...
    isequal(data.data.foo(1:2), data.columns{1}(1:2)) && isequal(data.data.bar, data.columns{2}) && ...
    isequal(data.data.Name_2, data.columns{3}) && isequal(data.data.empTy, data.columns{4}) && ...
    isnan(data.data.foo(3)), ...
    'csvload different delimiter, comment, fixed header, changed namefct');

f = fopen(fname, 'w');
fprintf(f, '1 2 3 4, x2\n');
fprintf(f, '2 3,  y \n');
fclose(f);
data = csvload(fname, 'trimws', false);
utexpect(isempty(data.header) && ...
    isequal(data.columns(:), {{[1 2 3 4] ; [2 3]} ; {' x2' ; '  y '}}) && ...
    isequal(data.coltypes, 'vs'), ...
    'csvload no trimming whitespace, no header, vector values');

quietdelete(fname);

end

function quietdelete(fname)

if exist(fname, 'file')
    delete(fname);
end

end