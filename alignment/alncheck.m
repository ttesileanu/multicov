function isaln = alncheck(alignment)
% ALNCHECK Return true if input is an alignment structure.
%   ALNCHECK(alignment) returns true if the argument has the correct
%   structure for an alignment.
%
%   See also ALNMAKE, BINCHECK, STATSCHECK.

% Tiberiu Tesileanu (2012-2014)

isaln = (genalncheck(alignment) && strcmp(alignment.type, 'character') && ...
    ismatrix(alignment.data) && ischar(alignment.data));

end