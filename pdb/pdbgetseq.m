function [seq, posnames] = pdbgetseq(pdb, pdbchain)
% PDBGETSEQ Get PDB sequence from PDB structure.
%   [seq, posnames] = PDBGETSEQ(pdbstruct, pdbchain) returns the PDB
%   sequence for the given chain of the given pdb structure (as loaded by
%   pdbread, for example). The function also returns the residue names, as
%   a cell array of strings.
%
% See also: PDBREAD.

% Tiberiu Tesileanu (2012, 2014)
%
% This function is loosely based on code by R. Ranganathan and K. Reynolds.

% some PDB structures don't have the 'Model' field and go directly to 'Atom'
try
    atom = pdb.Model.Atom;
catch ME
    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
        atom = pdb.Atom;
    else
        rethrow(ME);
    end
end

% mask to keep only the desired chain
chainmask = strcmpi({atom.chainID}, pdbchain);
atommasked = atom(chainmask);
% residue numbers
posnums = [atommasked.resSeq];
% take into account iCodes, if they exist
if isfield(atom, 'iCode')
    % in this case, combine residue numbers with iCodes
    icodes = {atommasked.iCode};
    keys = arrayfun(@(num, icode) ([int2str(num) icode{:}]), posnums, icodes, 'uniform', false);
else
    keys = arrayfun(@int2str, posnums, 'uniform', false);
end
% each residue has several atoms, so each number will appear several times
% we are assuming that names are unique
% we don't use the 'stable' flag because that only appeared in R2012a
[~, indsort0] = unique(keys, 'first');
indsort = sort(indsort0);
posnames = keys(indsort);
% [posnames, indsort] = unique(keys, 'stable');

% get the residues
residues = {atommasked(indsort).resName};

% transform amino acid residues to single-letter form
seq = blanks(length(residues));
for i = 1:length(residues)
    residue = residues{i};
    if length(residue) == 3
        seq(i) = aminolookup(residue);
    elseif length(residue) == 1
        seq(i) = residue;
    else
        error([mfilename ':badres'], ['Unrecognized residue "' residue '".']);
    end
end

end