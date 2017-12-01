function sf = splitstring(s, width)
% SPLITSTRING Format a string into a matrix of a given width.
%   sf = SPLITSTRING(s) returns a character matrix where the string s is
%   spread across several lines of width 78. The string is split only at
%   whitespace, if possible.
%
%   sf = SPLITSTRING(s, width) uses width 'width' instead.

% Tiberiu Tesileanu (2013-2014)

if nargin < 2
    width = 78;
end

sf = blanks(width);

crtpos = 1;
crtline = 1;
while crtpos <= length(s)
    % skip initial blanks
    substr = s(crtpos:end);
    sppos = find(~isspace(substr), 1, 'first');
    if isempty(sppos)
        return;
    else
        crtpos = crtpos + sppos - 1;
    end

    endpos = crtpos + width - 1;
    if endpos >= length(s)
        tocopy = length(s) - crtpos + 1;
    else
        substr = s(crtpos:endpos);
        sppos = find(isspace(substr), 1, 'last');
        if isempty(sppos) 
            tocopy = endpos - crtpos + 1;
        else
            tocopy = sppos - 1;
        end
    end
    if tocopy > 0
        sf(crtline, :) = blanks(width);
        sf(crtline, 1:tocopy) = s(crtpos:crtpos + tocopy - 1);
        crtline = crtline + 1;
    end
    crtpos = crtpos + tocopy;
end

end