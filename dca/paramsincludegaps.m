function [gapparams, changed] = paramsincludegaps(params, structure)
% PARAMSINCLUDEGAPS Include gap entries in the maximum entropy model
% couplings and fields.
%   gapparams = PARAMSINCLUDEGAPS(params) includes entries for the gaps in
%   the couplings from the maximum entropy model structure 'params'. These
%   are assumed to be in the 'gaps = 0' gauge, so all the inserted entries
%   are zero. This also updates the 'alphabets' field of the structure to
%   keep track of the change.
%
%   gapcouplings = PARAMSINCLUDEGAPS(couplings, structure) performs the
%   operation to the given coupling matrix using the given structure.
%
%   [..., changed] = PARAMSINCLUDEGAPS(...) returns 'changed' equal to true
%   if a change needed to be made, false otherwise.
%
% See also: INCLUDEGAPS, GETMAXENTPARAMS.

% Tiberiu Tesileanu (2014)

if ~paramscheck(params)
    if ~ismatrix(params) || ~isnumeric(params) || size(params, 1) ~= size(params, 2)
        error([mfilename ':badcouplings'], 'The first argument should either be a parameters structure or a square coupling matrix.');
    end
    if nargin < 2
        error([mfilename ':nostruct'], 'When the first argument is a coupling matrix, an alignment-like structure is needed as a second argument.');
    end
    [gapstructure, changed] = alphaincludegaps(structure);
    if changed
        binmap = cellfun(@(v) v(:), getbinmap(structure), 'uniform', false);
        gapbinmap = cellfun(@(v) v(:), getbinmap(gapstructure), 'uniform', false);
        gapparams = zeros(gapbinmap{end}(end));
        
        matchmap = cellfun(@(idxs, gapidxs) ...
            gapidxs(1+(length(idxs)~=length(gapidxs)):end), binmap, gapbinmap, ...
            'uniform', false);
        binvec = cell2mat(binmap(:));
        matchvec = cell2mat(matchmap(:));
        gapparams(matchvec, matchvec) = params(binvec, binvec);
    else
        gapparams = params;
    end
else
    [gapparams, changed] = alphaincludegaps(params);
    if changed
        gapparams.couplings = paramsincludegaps(params.couplings, params);
    else
        gapparams = params;
    end
end

end