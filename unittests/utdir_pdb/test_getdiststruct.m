function test_getdiststruct
% TEST_GETDISTSTRUCT Test getdiststruct.m.

% test this with mock coordinates
mock1 = {...
    struct('pos', [[-1, 0, 0] ; [0, 1, 0] ; [0, 0, 1]], 'type', {{'CA', 'CB', 'N'}}), ...   % com [-1/3, 1/3, 1/3]
    struct('pos', [[1, 0, 0] ; [0, 2, 0] ; [0, 0, 0]], 'type', {{'CA', 'CB', 'N'}}) ...     % com [1/3, 2/3, 0]
    struct('pos', [[1/2, 1/2, 0] ; [0, 1/2, 0] ; [0, 1/2, 1/2] ; [1/2, 1/2, 1/2]], ...
        'type', {{'CA', 'CB', 'N', 'S'}}) ...                                               % com [1/4, 1/2, 1/4]
};
mock2 = {...
    struct('pos', [[-1/3, 0, 1/3] ; [0, 1/3, 0] ; [0, 0, 1/3]], 'type', {{'CA', 'CB', 'N'}}), ...   % com [-1/9, 1/9, 2/9]
    struct('pos', [[0, 1/5, 0] ; [2/5, 2/3, 0] ; [0, 0, 1/4]], 'type', {{'CA', 'CB', 'N'}}) ...     % com [2/15, 13/45, 1/12]
};
dist = getdiststruct(mock1, mock2, 'method', 'ca');

utexpect(all(isfield(dist, {'distmat', 'validmask', 'refseq'})) && ...
    ismatrix(dist.distmat) && isnumeric(dist.distmat) && ...
    ismatrix(dist.validmask) && (isnumeric(dist.validmask) || islogical(dist.validsmask)) && ...
    isstruct(dist.refseq), ...
    'getdiststruct mock coordinates, output has proper fields');

utexpect(isvector(dist.refseq) && length(dist.refseq) == 2 && ...
    all(isfield(dist.refseq, {'seqdb', 'seqid', 'map'})), ...
    'getdiststruct refseq has proper size and fields');

utexpect(all(all(dist.validmask)), 'getdiststruct correct validmask');

distmat_exp = sqrt(...
   [  0     4   10/ 4   5/ 9  26/ 25 ; ...
      4     0    1/2   17/ 9  26/ 25 ; ...
    10/ 4  1/ 2   0    38/36  34/100 ; ...
     5/ 9 17/ 9 38/36    0    59/225 ; ...
    26/25 26/25 34/100 59/225   0]);
utexpect(all(abs(dist.distmat(:) - distmat_exp(:)) < 10*eps), ...
    'getdiststruct distance matrix from mock coordinates, ''method'' set to ''ca''');

dist_a = getdiststruct(mock1, mock2);
dist_a_exp = getdistmat([mock1 mock2]);
utexpect(all(abs(dist_a.distmat(:) - dist_a_exp(:)) < eps), ...
    'getdiststruct mock coordinates, default method');

% check that it works with PDBs
global multicov_ut_dir;

timer = tic;
pdb = pdbread(fullfile(multicov_ut_dir, 'data', '3TGI.pdb'));

dist1 = getdiststruct({pdb, 'E'}, 'method', 'ca');
coords_e = getpdbcoords(pdb, 'E');
distmat1_exp = getdistmat(coords_e, 'method', 'ca');
utexpect(all(all(dist.validmask)) && all(abs(dist1.distmat(:) - distmat1_exp(:)) < eps), ...
    'getdiststruct with single PDB chain, ''method'' ''ca''', timer);

timer = tic;
dist2 = getdiststruct({pdb, 'IE'}, 'method', 'ca', ...
    'refseq', struct('seqdb', {'uniprot', 'fake'}, ...
    'seqid', {'P1', 'faker'}, 'map', {[], 1:2:(2*length(coords_e)-1)}));
lens = arrayfun(@(r) length(r.map), dist2.refseq);
coords_i = getpdbcoords(pdb, 'I');
dist_ei_exp = zeros(length(coords_e), length(coords_i));
for i = 1:length(coords_e)
    posi = coords_e{i}.pos(strcmp(coords_e{i}.type, 'CA'), :);
    for j = 1:length(coords_i)
        posj = coords_i{j}.pos(strcmp(coords_i{j}.type, 'CA'), :);
        dist_ei_exp(i, j) = norm(posi - posj);
    end
end
echain_idxs = (lens(1)+1):sum(lens);
ichain_idxs = 1:lens(1);
utexpect(all(abs(flatten(dist2.distmat(echain_idxs, ichain_idxs) - dist_ei_exp)) < eps), ...
    'getdiststruct with multiple PDB chains', timer);

timer = tic;
distmat2_exp = distmat1_exp;
utexpect(all(abs(flatten(dist2.distmat(echain_idxs, echain_idxs) - distmat2_exp)) < eps), ...
    'getdiststruct within chain distances check', timer);

timer = tic;
dist2b = getdiststruct(pdb, 'method', 'ca', ...
    'refseq', struct('seqdb', {'fake', 'uniprot'}, ...
    'seqid', {'faker', 'P1'}, 'map', {1:2:(2*lens(2)-1), []}));
% this should be the same as dist2, but with the I and E chains reversed
reordering = [(lens(1)+1):sum(lens) 1:lens(1)];
dist2b.distmat(reordering, reordering) = dist2b.distmat;
dist2b.refseq([2 1]) = dist2b.refseq;
utexpect(isequal(dist2, dist2b), 'getdiststruct ordering of PDB chains', timer);

% check mixed PDB / coords input
timer = tic;
dist2c = getdiststruct({pdb, 'I'}, coords_e, 'method', 'ca');
utexpect(all(abs(dist2.distmat(:) - dist2c.distmat(:)) < eps), ...
    'getdiststruct mixed PDB / coords input', timer);

% check refseq override
utexpect(isequal(dist2.refseq(1), dist2c.refseq(1)) && ...
    strcmp(dist2.refseq(2).seqdb, 'fake') && ...
    strcmp(dist2.refseq(2).seqid, 'faker') && ...
    isequal(dist2.refseq(2).map(:), flatten(1:2:(2*lens(2)-1))), ...
    'getdiststruct partially override refseq');

% check taking refseq from PDB file
timer = tic;
pdb2 = pdbread(fullfile(multicov_ut_dir, 'data', '4TNA.pdb'));

dist3 = getdiststruct({pdb2, 'A'});
[~, posnames2] = pdbgetseq(pdb2, 'A');
utexpect(isequal(dist3.refseq.map(:), posnames2(:)), ...
    'getdiststruct with refseq mapping from PDB residue names', timer);

end