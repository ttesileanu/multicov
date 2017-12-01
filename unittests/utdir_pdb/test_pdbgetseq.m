function test_pdbgetseq
% TEST_PDBGETSEQ Test pdbgetseq.m

global multicov_ut_dir;

timer = tic;
pdb = pdbread(fullfile(multicov_ut_dir, 'data', '3TGI.pdb'));
pdb_chain = 'E';

[seq, names] = pdbgetseq(pdb, pdb_chain);
seq0 = pdb.Sequence(strcmp({pdb.Sequence.ChainID}, pdb_chain)).Sequence;
chainmask = strcmp({pdb.Model.Atom.chainID}, pdb_chain);
nums0 = arrayfun(@int2str, [pdb.Model.Atom(chainmask).resSeq], 'uniform', false);
icodes0 = {pdb.Model.Atom(chainmask).iCode};
names0unsorted = cellfun(@(n, i) [n i], nums0, icodes0, 'uniform', false);
[~, tmp] = unique(names0unsorted, 'first');
names0 = names0unsorted(sort(tmp));
utexpect(strcmp(seq, seq0) && isequal(names, names0), 'pdbgetseq', timer);

end