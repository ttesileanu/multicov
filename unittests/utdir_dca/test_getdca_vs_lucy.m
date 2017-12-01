function test_getdca_vs_lucy
% TEST_GETDCA_VS_LUCY Test that getdca.m returns results identical to
% Lucy's code.

global multicov_ut_dir;

% load the alignment and Lucy's results
orig = open(fullfile(multicov_ut_dir, 'data', 'lucy_dca_pdz_small.mat'));

% get our own DCA matrix
timer = tic;
dca_noseqw = getdca(orig.alignment);
utexpect(all(abs(dca_noseqw.cmat(:) - orig.DCAmat_noseqw(:)) < eps), ...
    'getdca no sequence weights', timer);

timer = tic;
seqw = estimateseqw(orig.alignment, 0.7);
utexpect(all(abs(seqw(:) - orig.seqw(:)) < eps), ...
    'getdca estimateseqw result', timer);

timer = tic;
alignment = orig.alignment;
alignment.seqw = seqw;
dca = getdca(alignment);
utexpect(all(abs(dca.cmat(:) - orig.DCAmat(:)) < eps), ...
    'getdca with sequence weights', timer);

end