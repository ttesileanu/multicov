function [type, value] = guesstype(s)
% GUESSTYPE Guess the type of a string.
%   [type, value] = GUESSTYPE(s) tries to identifying the type of a string.
%   Here are the recognized types, in the order in which they are detected:
%     'x':  empty, returned if s == ''; value returned is also ''.
%     'l':  logical, returned if strtrim(lower(s)) is 'true' or 'false';
%           value is logical true or false
%     'n':  number; value is the numeric equivalent of s
%     'v':  space-separated vector of numbers; value is the vector
%     's':  if all else fails; value is s

% Tiberiu Tesileanu (2013-2014)

if ~ischar(s)
    error([mfilename ':badarg'], 'The argument should be a string.');
end

value = s;
if isempty(s)
    type = 'x';
    return;
end

% check for logical values
s1 = strtrim(s);
if strcmpi(s1, 'true')
    type = 'l';
    value = true;
    return;
end
if strcmpi(s1, 'false')
    type = 'l';
    value = false;
    return;
end

% try to convert to a number
n = str2double(s);
if isnan(n)
    % this either wasn't a number to start with, or it was an explicit nan
    ok = strcmpi(s1, 'nan');
else
	ok = true;
end

if ok
    value = n;
    type = 'n';
    return;
end

% it wasn't a number... try a vector
% split the string up at white space
sppos = [0 find(isstrprop(s, 'wspace')) length(s)+1];
elems = arrayfun(@(i) s((sppos(i) + 1):(sppos(i+1)-1)), 1:(length(sppos)-1), 'uniform', false);
% get rid of empty elements, that may appear because of repeated white space
elems = elems(~cellfun(@isempty, elems));
% check that all elements are convertible to numbers
nelems = cellfun(@str2double, elems);
nanmask = strcmpi(elems, 'nan');
ok = ~any(isnan(nelems(~nanmask)));

if ok
    value = nelems;
    type = 'v';
    return;
end

% all else failed -- it's just a string
type = 's';

end