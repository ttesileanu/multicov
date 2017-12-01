function test_dca_vs_lucy
% TEST_DCA_VS_LUCY Test that the DI and spectral norm reductions of the DCA
% matrix match those from Lucy's code.

tol = 1e-12;

global multicov_ut_dir;

% load the alignment and Lucy's results
orig = open(fullfile(multicov_ut_dir, 'data', 'lucy_dca_pdz_small.mat'));

% get DCA with and without sequence weights
timer = tic;
dca_noseqw = getdca(orig.alignment);
seqw = estimateseqw(orig.alignment, 0.7);
alignment = orig.alignment;
alignment.seqw = seqw;
dca = getdca(alignment);

params_noseqw = getmaxentparams(dca_noseqw);
params = getmaxentparams(dca);

% calculate DI
dimat_noseqw = reducecouplings(params_noseqw, 'di', 'freq1', dca_noseqw.freq1);
dimat = reducecouplings(params, 'di', 'freq1', dca.freq1);
utexpect(all(abs(dimat_noseqw(:) - flatten(zeroblockdiag(orig.DImat_noseqw, 1))) < tol) && ...
    all(abs(dimat(:) - flatten(zeroblockdiag(orig.DImat, 1))) < tol), ...
    'reducecouplings DI calculation vs. Lucy, with/without sequence weights', timer);

% calculate spectral norm
timer = tic;
specmat_noseqw = reducecouplings(params_noseqw, 'norm', 'norm', 'spectral');
specmat = reducecouplings(params, 'norm', 'norm', 'spectral');
utexpect(all(abs(specmat_noseqw(:) - flatten(zeroblockdiag(orig.Specmat_noseqw, 1))) < tol) && ...
    all(abs(specmat(:) - flatten(zeroblockdiag(orig.Specmat, 1))) < tol), ...
    'reducecouplings spectral norm vs. Lucy, with/without sequence weights', timer);

% XXX write code for and test MI

end