function test_getpdbcoords
% TEST_GETPDBCOORDS Test getpdbcoords.m.

global multicov_ut_dir;

timer = tic;
pdb = pdbread(fullfile(multicov_ut_dir, 'data', '3TGI.pdb'));
chain = 'E';
coords1 = getpdbcoords(pdb, chain);
chainmask = find(strcmp({pdb.Model.Atom.chainID}, chain));
resnames = arrayfun(@int2str, [pdb.Model.Atom(chainmask).resSeq], 'uniform', false);
icodes = {pdb.Model.Atom(chainmask).iCode};
posnames0 = cellfun(@(s1, s2) [s1 s2], resnames, icodes, 'uniform', false);
% can't use 'stable' before R2012a, better to avoid it
% posnames = unique(posnames0, 'stable');
[~, tmp] = unique(posnames0, 'first');
posnames = posnames0(sort(tmp));
n = length(posnames);
coords_exp = cell(n, 1);
for i = 1:n
    resmask = strcmp(posnames0, posnames{i});
    thismask = chainmask(resmask);
    coords_exp{i}.pos = [pdb.Model.Atom(thismask).X ; pdb.Model.Atom(thismask).Y ; pdb.Model.Atom(thismask).Z]';
    coords_exp{i}.type = {pdb.Model.Atom(thismask).AtomName}';
end
utexpect(isequal(coords1, coords_exp), 'getpdbcoords', timer);

end