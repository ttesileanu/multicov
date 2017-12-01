function test_getsca_vs_rama
% TEST_GETSCA_VS_RAMA Test that getsca.m behaves as Rama Ranganathan's SCA
% in cases where they are both applicable. This tests versions 4.0, 5.0,
% and the binary approximation of SCA. Also tests that getconservation.m
% gives the same results as Rama Ranganathan's cons.m when used with the KL
% divergence.

global multicov_ut_dir;

% load the alignment and Rama's results
orig = open(fullfile(multicov_ut_dir, 'data', 'rama_sca_pdz_small.mat'));

% XXX we are assuming differences are fine if they're within 1e-10

tol = 1e-10;

timer = tic;
sca4 = getsca(orig.alignment);
utexpect(all(abs(sca4.cmat(:) - orig.Csca4(:)) < tol), ...
    'getsca vs. SCA 4.0', timer);

timer = tic;
sca5 = getsca(orig.alignment, 'method', 'projection');
utexpect(all(abs(sca5.cmat(:) - orig.Csca5(:)) < tol), ...
    'getsca vs. SCA 5.0', timer);

timer = tic;
sca_bin = getsca(orig.alignment, 'method', 'binary');
utexpect(all(abs(sca_bin.cmat(:) - orig.Csca_bin(:)) < tol), ...
    'getsca vs. binary SCA', timer);

% we also test the conservation here... not quite the proper place, but
% close enough
cons = getconservation(orig.alignment, 'method', 'kl', 'gaps', true);
utexpect(all(abs(cons(:) - orig.consvec(:)) < tol), ...
    'getsca getconservation vs. SCA code');

end