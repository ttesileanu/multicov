function test_csvsave
% TEST_CSVSAVE Test csvsave.m.

global multicov_ut_dir;

fname = fullfile(multicov_ut_dir, 'temp', 'csvsave_test.csv');

csvsave(fname, [2 3 4 ; 5 6 7], 'header', {'a', 'b', 'c'});
f = fopen(fname);
utexpect(f ~= -1, 'csvsave test file existence');
cdata = textscan(f, '%s', 'delimiter', '\n');
cdata = cdata{1};
fclose(f);
utexpect(isequal(cdata, {'a,b,c' ; '2,3,4' ; '5,6,7'}), ...
    'csvsave numeric matrix with header');

csvsave(fname, {'foo', 'bar' ; 'aleph', 'beth' ; 'john', 'mary'});
f = fopen(fname);
cdata = textscan(f, '%s', 'delimiter', '\n');
cdata = cdata{1};
fclose(f);
utexpect(isequal(cdata, {'foo,bar' ; 'aleph,beth' ; 'john,mary'}), ...
    'csvsave cell matrix');

csvsave(fname, {[1 2 3 4 5], [2 4 6 8 10]}, 'header', {'x', '2x'}, 'delimiter', ' ');
f = fopen(fname);
cdata = textscan(f, '%s', 'delimiter', '\n');
cdata = cdata{1};
fclose(f);
utexpect(isequal(cdata, {'x 2x' ; '1 2' ; '2 4' ; '3 6' ; '4 8' ; '5 10'}), ...
    'csvsave cell of columns, space delimiter');

quietdelete(fname);

end

function quietdelete(fname)

if exist(fname, 'file')
    delete(fname);
end

end