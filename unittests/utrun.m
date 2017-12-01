function [result, stats] = utrun(varargin)
% UTRUN Run the unit tests in the current folder and its subfolders.
%   UTRUN runs all the unit tests found in the current folder and in its
%   subfolders. Each test must be a function with a name starting with the
%   characters 'test_'. The tests must use the utexpect function to notify
%   the system of success or failure. Any uncaught exception thrown by the
%   test functions is interpreted as a failure.
%
%   [result, stats] = UTRUN returns result = true or false, depending on
%   whether all the tests passed or not, and a stats structure containing
%   statistics of the run:
%    'ncases':  Total number of test cases.
%    'nchecks': Total number of times utexpect was run.
%    'fcases':  Fraction of successful cases.
%    'fchecks': Fraction of successful checks.
%
%   This function can take the following options:
%    'cpppolicy' <s>
%       How to handle tests of C++ extensions when these aren't found.
%         'fail':   fail test if a C++ extension can't be found
%         'warn':   issue a warning but pass the test
%         'pass':   silently pass the test
%       (default: 'warn')
%    'dirchoice' <s>
%       Regex describing the unit test folders to be run.
%       (default: all)
%    'verbosity' <v>
%       Define the verbosity level.
%         -2:   no screen output
%         -1:   only output at end of all test cases
%          0:   only output at end of a test case
%          1:   output at each check, but only on failure
%          2:   output at each check
%       (default: 1)
%
% See also: UTEXPECT.

% Tiberiu Tesileanu (2014)

old_warning_state = warning;
warning('on'); %#ok<WNON>

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

% parameters
parser.addParamValue('verbosity', 1, @(x) x == floor(x));
parser.addParamValue('dirchoice', '.*', @(s) isempty(s) || (ischar(s) && isvector(s)));
parser.addParamValue('cpppolicy', 'warn', @(s) ismember(s, {'fail', 'warn', 'pass'}));

% parse
parser.parse(varargin{:});
params = parser.Results;

% set up the global variables updated by utexpect
global multicov_ut_verbosity;
global multicov_ut_test_success;
global multicov_ut_nchecks;
global multicov_ut_ncheck_successes;
global multicov_ut_dir;
global multicov_ut_cpppolicy;

multicov_ut_verbosity = params.verbosity;
multicov_ut_cpppolicy = params.cpppolicy;
multicov_ut_nchecks = 0;
multicov_ut_ncheck_successes = 0;

% keep track of the current directory to return to it at the end
PWD = pwd;

% move to the directory in which utrun is
multicov_ut_dir = fileparts(mfilename('fullpath'));
cd(multicov_ut_dir);

% check that we ran setup_paths, and get all the files in this directory
% and its subdirectories
try
    files = findfiles;
catch me
    if strcmp(me.identifier, 'MATLAB:UndefinedFunction')
        error('utrun:setuperror', 'Function ''findfiles'' not defined -- did you run setup_paths?');
    else
        rethrow(me);
    end
end

% identify the test-case files
get_extension = @(s) getnth(3, @fileparts, s);
is_test = @(s) length(s) > 5 && strcmp(s(1:5), 'test_');
testmask = cellfun(@(s) is_test(getnth(2, @fileparts, s)) && strcmp(get_extension(s), '.m'), files);
test_files = files(testmask);
paths = cellfun(@fileparts, test_files, 'uniform', false);
[~, tests] = cellfun(@fileparts, test_files, 'uniform', false);

ncase_successes = 0;

timer = tic;

% run the tests
ntotaltests = 0;
for i = 1:length(tests)
    if isempty(regexp(paths{i}, params.dirchoice, 'once'))
        if params.verbosity >= 2
            % skip this test
            disp(['UT case:  ' test_files{i} ': skipped because of dirchoice flag.']);
        end
        continue;
    end
    
    ntotaltests = ntotaltests + 1;
    
    multicov_ut_test_success = true;
    % go to the right directory and run the test
    test_timer = tic;
    cd(paths{i});
    try
        eval([tests{i} ';']);
    catch me
        % handle exceptions
        multicov_ut_test_success = false;
        if params.verbosity >= 1
            bt_state = warning('query', 'backtrace');
            warning('off', 'backtrace');
            warning('test:fail', ['Exception in ' test_files{i} ...
                '. Here are the details: ' me.identifier ': ' me.message]);
            for j = 1:length(me.stack)-1
                warning('test:fail', [' trace: ' me.stack(j).name ' (line ' ...
                    int2str(me.stack(j).line) ')']);
            end
            warning(bt_state);
        end
    end
    
    test_dt = toc(test_timer);
    % display the outcome of the test
    if params.verbosity >= 0
        msg = ['UT case:  ' test_files{i} ': '];
        if test_dt > 1
            timer_msg = [' (' num2str(test_dt) ' seconds)'];
        else
            timer_msg = '';
        end
        if multicov_ut_test_success
            msg = [msg 'pass' timer_msg]; %#ok<AGROW>
            disp(msg);
        else
            msg = [msg 'FAIL!' timer_msg]; %#ok<AGROW>
            bt_state = warning('query', 'backtrace');
            warning('off', 'backtrace');
            warning('test:fail', msg);
            warning(bt_state);
        end
    end
    
    % count successes
    if multicov_ut_test_success
        ncase_successes = ncase_successes + 1;
    end
    
    % go back to the base directory
    cd(multicov_ut_dir);
end

% output success stats to screen
if params.verbosity >= -1
    disp(' ');
    if ncase_successes == ntotaltests
        disp(['All ' int2str(ntotaltests) ' tests PASSED.']);
    else
        bt_state = warning('query', 'backtrace');
        warning('off', 'backtrace');
        warning('test:stats', 'Some FAILURES, %d of %d tests passed.', ncase_successes, ntotaltests);
        warning(bt_state);
    end
    
    disp(['Total time: ' num2str(toc(timer)) ' seconds.']);
end

% return some success stats
result = (ncase_successes == ntotaltests);

stats.nfiles = length(files);
stats.ncases = ntotaltests;
stats.nchecks = multicov_ut_nchecks;
stats.fcases = ncase_successes / stats.ncases;
stats.fchecks = multicov_ut_ncheck_successes / stats.nchecks;

% go back to the directory we were in when the function was called
cd(PWD);

warning(old_warning_state);

end