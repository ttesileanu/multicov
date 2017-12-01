function [res, changed] = excludegaps(orig)
% EXCLUDEGAPS Exclude gaps from consideration for frequency vectors,
% alignments, statistics structures, or maximum entropy model parameter
% structures.
%   res = EXCLUDEGAPS(orig) removes gaps from the argument, calling
%   alphaexcludegaps, alnexcludegaps, statsexcludegaps, or paramsexcludegaps
%   depending on its form.
%
%   [res, changed] = EXCLUDEGAPS(orig) returns 'changed' equal to true if a
%   change needed to be made, false otherwise.
%
% See also: ALPHAEXCLUDEGAPS, ALNEXCLUDEGAPS, STATSEXCLUDEGAPS,
% PARAMSEXCLUDEGAPS, INCLUDEGAPS.

% Tiberiu Tesileanu (2014)

if alncheck(orig) || bincheck(orig)
    [res, changed] = alnexcludegaps(orig);
elseif statscheck(orig)
    [res, changed] = statsexcludegaps(orig);
elseif paramscheck(orig)
    [res, changed] = paramsexcludegaps(orig);
elseif iscell(orig) || ischar(orig) || isstruct(orig)
    [res, changed] = alphaexcludegaps(orig);
else
    error([mfilename ':badarg'], 'This function works with character or binary alignments, statistics structures, maxent structures, alphabets, or cell arrays of alphabets.');
end

end