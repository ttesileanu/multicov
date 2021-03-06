function [newalphas, changed] = alphaincludegaps(alphas)
% ALPHAINCLUDEGAPS Turn a list of alphabets into gappy alphabets.
%   newalphas = ALPHAINCLUDEGAPS(alphas) takes in a single alphabet or a
%   cell array of alphabets, and returns it with 'gapped' alphabets instead
%   (adding the letters 'gap' before each of the names).
%
%   newalign = ALPHAINCLUDEGAPS(structure) does this to the 'alphabets'
%   field of a structure.
%
%   [newalign, changed] = ALPHAINCLUDEGAPS(...) also returns 'changed'
%   equal to true if a change was performed, false otherwise.
%
% See also: ALPHAGETLETTERS, INCLUDEGAPS, ALNINCLUDEGAPS.

% Tiberiu Tesileanu (2014)

if isstruct(alphas)
    if ~isfield(alphas, 'alphabets')
        error([mfilename ':badarg'], 'When the argument is a structure, it should have a field named ''alphabets''.');
    end
    newalphas = alphas;
    [newalphas.alphabets, changed] = alphaincludegaps(alphas.alphabets);
    return;
elseif ~iscell(alphas)
    nocell = true;
    alphas = {alphas};
else
    nocell = false;
end

newalphas = alphas;

% turn the alphabets into gapped alphabets
changed = false;
for a = 1:length(alphas)
    % do nothing if this is already a 'gapped' alphabeet
    if length(alphas{a}) < 3 || ~strcmp(alphas{a}(1:3), 'gap')        
        newalphas{a} = ['gap' alphas{a}];
        changed = true;
    end
end

if nocell
    newalphas = newalphas{1};
end

end