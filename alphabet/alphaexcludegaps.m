function [newalphas, changed] = alphaexcludegaps(alphas)
% ALPHAEXCLUDEGAPS Turn a list of alphabets into non-gappy alphabets.
%   newalphas = ALPHAEXCLUDEGAPS(alphas) takes in a single alphabet or a
%   cell array of alphabets, and returns it after replacing all 'gapped'
%   alphabets with their ungapped counterparts (i.e., removing the 'gap'
%   prefix from the alphabets that have it).
%
%   newalign = ALPHAEXCLUDEGAPS(structure) does this to the 'alphabets'
%   field of a structure.
%
%   [newalign, changed] = ALPHAEXCLUDEGAPS(...) returns 'changed' equal to
%   true if a change was performed, false otherwise.
%
% See also: ALPHAGETLETTERS, INCLUDEGAPS, EXCLUDEGAPS.

% Tiberiu Tesileanu (2014)

if isstruct(alphas)
    if ~isfield(alphas, 'alphabets')
        error([mfilename ':badarg'], 'When the argument is a structure, it should have a field named ''alphabets''.');
    end
    newalphas = alphas;
    [newalphas.alphabets, changed] = alphaexcludegaps(alphas.alphabets);
    return;
elseif ~iscell(alphas)
    nocell = true;
    alphas = {alphas};
else
    nocell = false;
end

newalphas = alphas;

% turn the alphabets into ungapped alphabets
changed = false;
for a = 1:length(alphas)
    % do nothing if this is already an 'ungapped' alphabeet
    if length(alphas{a}) > 3 && strcmp(alphas{a}(1:3), 'gap')        
        newalphas{a} = alphas{a}(4:end);
        changed = true;
    end
end

if nocell
    newalphas = newalphas{1};
end

end