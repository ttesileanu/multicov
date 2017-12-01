function res = tofieldname(s)
% TOFIELDNAME Turn the string into a proper structure field name.
%   res = TOFIELDNAME(s) turns the string into a proper structure field
%   name by making sure that it starts with a letter, and replacing any
%   invalid characters with underscores.
%
%   This is similar to Matlab's genvarname, but that function treats
%   symbols awkwardly (by transforming them to ASCII codes). This function
%   instead takes a more user-friendly approach.
%
% See also: GENVARNAME.

% Tiberiu Tesileanu (2013-2014)

tmp = find(isletter(s));
if isempty(tmp)
    error([mfilename ':bads'], ['The string "' s '" cannot be turned into a proper field name because it contains no letters.']);
end

% start from the first letter
res = s(tmp(1):end);

% turn everything that's not a letter, a digit, or an underscore into
% underscores
res(~isletter(res) & ~isstrprop(res, 'digit')) = '_';

end