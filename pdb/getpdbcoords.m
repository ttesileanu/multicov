function coords = getpdbcoords(pdb, chain)
% GETPDBCOORDS Get the coordinates of the atoms in the PDB structure.
%   coords = GETPDBCOORDS(pdb, chain) returns a cell array of coordinates
%   and types for the atoms in each of the residues of the PDB chain. Each
%   cell of coords is a structure with a field 'pos' containing a K x 3
%   matrix of positions (one row per atom), and a field 'type' containing a
%   cell array of strings giving the atom types.

% Tiberiu Tesileanu (2012, 2014)

% some PDB structures have directly pdb.Atom instead of pdb.Model.Atom
try
    atom = pdb.Model.Atom;
catch ME
    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
        atom = pdb.Atom;
    else
        throw(ME);
    end
end

% mask to keep only the correct chain
chainidxs = find(strcmp({atom.chainID}, chain));

coords = {};
lastres = '';
currentres = [];
currentname = {};
for i = 1:length(chainidxs)
    j = chainidxs(i);
    res = [num2str(atom(j).resSeq) atom(j).iCode];
    name = atom(j).AtomName;
    if ~strcmp(res, lastres)
        if i > 1
            % starting new residue
            coords = [coords ; struct('pos', currentres, 'type', {currentname})]; %#ok<AGROW>
            currentres = [];
            currentname = {};
        end
        lastres = res;
    end
    
    newatom = [atom(j).X atom(j).Y atom(j).Z];
    currentres = [currentres ; newatom]; %#ok<AGROW>
    currentname = [currentname ; name]; %#ok<AGROW>
end

% add the last residue
coords = [coords ; struct('pos', currentres, 'type', {currentname})];

end