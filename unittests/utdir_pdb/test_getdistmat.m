function test_getdistmat
% TEST_GETDISTMAT Test getdistmat.m.

mock = {...
    struct('pos', [[-1, 0, 0] ; [0, 1, 0] ; [0, 0, 1]], 'type', {{'CA', 'CB', 'N'}}), ...   % com [-1/3, 1/3, 1/3]
    struct('pos', [[1, 0, 0] ; [0, 2, 0] ; [0, 0, 0]], 'type', {{'CA', 'CB', 'N'}}) ...     % com [1/3, 2/3, 0]
};

dist1 = getdistmat(mock, 'method', 'min');
dist_min = [0 1 ; 1 0];
utexpect(abs(dist1(:) - dist_min(:)) < eps, 'getdistmat with coords, min');

dist2 = getdistmat(mock, 'method', 'gcm');
dist_gcm_12 = norm([-1/3, 2/3, 0] - [-1/3, 1/3, 1/3]);
dist_gcm = [0 dist_gcm_12 ; dist_gcm_12 0];
utexpect(abs(dist2(:) - dist_gcm(:)) < eps, 'getdistmat with coords, gcm');

global multicov_ut_dir;

timer = tic;
pdb = pdbread(fullfile(multicov_ut_dir, 'data', '3TGI.pdb'));
chain = 'E';
dist3 = getdistmat({pdb, chain}, 'method', 'ca');
chainmask = strcmp({pdb.Model.Atom.chainID}, chain);
camask = (chainmask & strcmp({pdb.Model.Atom.AtomName}, 'CA'));
capos = [pdb.Model.Atom(camask).X ; pdb.Model.Atom(camask).Y ; pdb.Model.Atom(camask).Z]';
dist_ca = squareform(pdist(capos));
utexpect(abs(dist3(:) - dist_ca(:)) < eps, 'getdistmat with pdb, ca', timer);

end