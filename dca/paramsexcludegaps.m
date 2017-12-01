function [ungapparams, changed] = paramsexcludegaps(params, structure)
% PARAMSEXCLUDEGAPS Exclude gap entries from the maximum entropy model
% couplings and fields.
%   ungapparams = PARAMSEXCLUDEGAPS(params) excludes entries for the gaps
%   from the couplings in the maximum entropy model structure 'params'.
%   This also updates the 'alphabets' field of the structure to keep track
%   of the change.
%
%   ungapcouplings = PARAMSEXCLUDEGAPS(couplings, structure) performs the
%   operation to the given coupling matrix using the given structure.
%
%   [..., changed] = PARAMSEXCLUDEGAPS(...) returns 'changed' equal to true
%   if a change needed to be made, false otherwise.
%
% See also: EXCLUDEGAPS, GETMAXENTPARAMS.

% Tiberiu Tesileanu (2014)

if ~paramscheck(params)
    if ~ismatrix(params) || ~isnumeric(params) || size(params, 1) ~= size(params, 2)
        error([mfilename ':badcouplings'], 'The first argument should either be a parameters structure or a square coupling matrix.');
    end
    if nargin < 2
        error([mfilename ':nostruct'], 'When the first argument is a coupling matrix, an alignment-like structure is needed as a second argument.');
    end
    [ungapstructure, changed] = alphaexcludegaps(structure);
    if changed
        binmap = cellfun(@(v) v(:), getbinmap(structure), 'uniform', false);
        ungapbinmap = cellfun(@(v) v(:), getbinmap(ungapstructure), 'uniform', false);
        ungapparams = zeros(ungapbinmap{end}(end));
        matchmap = cellfun(@(idxs, gapidxs) ...
            gapidxs(1+(length(idxs)~=length(gapidxs)):end), ungapbinmap, binmap, ...
            'uniform', false);
        ungapbinvec = cell2mat(ungapbinmap(:));
        matchvec = cell2mat(matchmap(:));
        ungapparams(ungapbinvec, ungapbinvec) = params(matchvec, matchvec);
    else
        ungapparams = params;
    end
else
    [ungapparams, changed] = alphaexcludegaps(params);
    if changed
        ungapparams.couplings = paramsexcludegaps(params.couplings, params);
    else
        ungapparams = params;
    end
end

end