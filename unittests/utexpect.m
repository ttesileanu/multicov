function utexpect(condition, msg_prefix, timer)
% UTEXPECT Output success/fail message based on condition.
%   UTEXPECT(condition, msg_prefix) outputs an appropriate message
%   depending on whether the condition is true or false.

% Tiberiu Tesileanu (2014)

global multicov_ut_verbosity;
global multicov_ut_test_success;
global multicov_ut_nchecks;
global multicov_ut_ncheck_successes;

multicov_ut_nchecks = multicov_ut_nchecks + 1;
if condition
    multicov_ut_ncheck_successes = multicov_ut_ncheck_successes + 1;
end

verb = multicov_ut_verbosity;

% check if we were asked to time the task
timer_text = '';
if nargin > 2
    dt = toc(timer);
    % but don't display this for very short times
    if dt >= 0.1
        timer_text = ['   (calculation took ' num2str(dt) ' seconds)'];
    end
end
msg0 = ['  UT check: ' msg_prefix ': '];
if ~condition
    multicov_ut_test_success = false;
    if verb >= 1
        msg = [msg0 'FAIL!' timer_text];
        % display in red for emphasis
        % but avoid displaying the backtrace
        bt = warning('query', 'backtrace');
        warning('off', 'backtrace');
        warning('expect:fail', msg);
        % reset backtrace display to whatever settings there were before
        warning(bt);
    end
elseif verb >= 2
    msg = [msg0 'pass' timer_text];
    disp(msg);
end

end