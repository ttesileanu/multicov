function [res, changed] = includegaps(orig)
% INCLUDEGAPS Add gaps to frequency vectors, alignments, statistics
% structures, or maximum entropy model parameter structures.
%   res = INCLUDEGAPS(orig) adds gaps to the argument, calling
%   alphaincludegaps, alnincludegaps, statsincludegaps, or paramsincludegaps
%   depending on its form.
%
%   [res, changed] = INCLUDEGAPS(orig) returns 'changed' equal to true if a
%   change needed to be made, false otherwise.
%
% See also: ALPHAINCLUDEGAPS, ALNINCLUDEGAPS, STATSINCLUDEGAPS,
% PARAMSINCLUDEGAPS, EXCLUDEGAPS.

% Tiberiu Tesileanu (2014)

if alncheck(orig) || bincheck(orig)
    [res, changed] = alnincludegaps(orig);
elseif statscheck(orig)
    [res, changed] = statsincludegaps(orig);
elseif paramscheck(orig)
    [res, changed] = paramsincludegaps(orig);
elseif iscell(orig) || ischar(orig) || isstruct(orig)
    [res, changed] = alphaincludegaps(orig);
else
    error([mfilename ':badarg'], 'This function works with character or binary alignments, statistics structures, maxent structures, alphabets, or cell arrays of alphabets.');
end

end