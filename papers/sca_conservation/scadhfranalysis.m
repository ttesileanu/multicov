% This script performs the analysis of the SCA method on the DHFR
% alignments, as described in the paper. This is also the script used to
% generate the graphs.
%
% The basic purpose of this script is to show how the sector found in
% Reynolds et al. 2011 compares to the top conserved positions, how
% the comparison changes if we change the alignment or the threshold, etc.
%
% There are some options that can be set by defining some variables before
% running the first section of the script. Leave these variables undefined
% to get the default behavior.
%   'aln_choice'
%       Choice of alignment to use. Can be
%         'hhblits_e3': HHblits alignment, e-value: e3
%         'hhblits_e10': HHblits alignment, e-value: e10
%         'rama': the original alignment from Reynolds et al. 2011,
%         'pfam': raw Pfam alignment, or
%       (default: 'hhblits_e10')
%   'cons_bkg_choice'
%       Choice of background distribution to be used for conservation
%       calculations (this only affects the KL conservation). Options are
%         'rama':    the usual distribution used by SCA (use 'getdbbkg' to
%                    see what this is)
%         'uniform': uniform distribution
%       or it can be a vector of 20 numbers that add up to 1 (the
%       frequencies of the amino acids in the order returned by
%       alphagetletters('protein', 'nogap')
%       (default: 'rama')
%   'cutoffrad'
%       Cutoff radius (in angstroms) to use for finding which residues
%       'touch' which surface residues. As in Reynolds et al., 'touching'
%       is defined as any atom being inside one of the spheres centered on
%       the four peptide bond atoms where the LOV2 insertion takes place,
%       where the radius of the spheres is given by 'cutoffrad'.
%       (default: 4)
%   'filtering'
%       Whether to filter the alignment (cap gaps at 20% per row and 20% per
%       column, and filter out sequences that are more than 80% identical).
%       (default: false)
%   'make_figures'
%       Whether to make the figures.
%       (default: true)
%   'sca_bkg_choice'
%       A choice of background distribution to be used with SCA. See
%       'cons_bkg_choice' for choices.
%       (default: 'rama')
%   'sca_type'
%       Type of SCA matrix to use. Options:
%        'reduction':   reduction method SCA. (version 4.0)
%        'projection':  projection method SCA. (version 5.0)
%        'binary':      binary approximation SCA.
%        'reynolds':    the reduction method used in Reynolds et al. (SCA 4.5)
%       (default: 'projection')
%   'sec_fraction'
%       Only when 'uxse_top_ev' is true: fraction of protein to be used for
%       the sector. The protein size is assumed to be equal to the width of
%       the alignment. If sec_fraction is set to a number larger than or
%       equal to 1, it is assumed to be equal to the number of residues
%       that should be taken up by the sector.
%       (default: 0.25)
%   'sub_rnd'
%       Whether to subtract from the SCA matrix for the alignment the
%       average SCA matrix obtained for randomized versions of the
%       alignment. Reynolds et al. uses this.
%       (default: false)
%   'use_top_ev'
%       When this is true, only the top eigenvector of the SCA matrix is
%       used. When it is false, the top 5 eigenvectors are used as in
%       Renyolds et al.; in this latter case, the Reynolds et al. procedure
%       and thresholds are used, and 'sec_fraction' is ignored.
%       (default: true)
%   'zlimit'
%       Threshold to use for the experimental z-score.
%       (default: 2)

% Tiberiu Tesileanu (2013-2014)

%% I. Information we get from the alignment alone

%% 1. Reproducing the results in Reynolds et al.

% handle some options
if ~exist('aln_choice', 'var')
    aln_choice = 'hhblits_e10';
end
if ~exist('make_figures', 'var')
    make_figures = true;
end
if ~exist('filtering', 'var')
    filtering = false;
end
if ~exist('zlimit', 'var')
    zlimit = 2;
end
if ~exist('sub_rnd', 'var')
    sub_rnd = false;
end
if ~exist('cutoffrad', 'var')
    cutoffrad = 4;
end
if ~exist('sec_fraction', 'var')
    sec_fraction = 0.25;
end
if ~exist('sca_type', 'var')
    sca_type = 'projection';
end
if ~exist('use_top_ev', 'var')
    use_top_ev = true;
end

if ~exist('cons_bkg_choice', 'var')
    cons_bkg_choice = 'rama';
else
    if strcmp(cons_bkg_choice, 'uniform')
        cons_bkg_choice = ones(20, 1)/20;
    elseif ~strcmp(sca_bkg_choice, 'rama') && ~isnumeric(sca_bkg_choice)
        error('script:consbkg', 'cons_bkg_choice should be ''rama'', ''uniform'', or a vector of 20 frequencies');
    end
end
if ~exist('sca_bkg_choice', 'var')
    sca_bkg_choice = 'rama';
else
    if strcmp(sca_bkg_choice, 'uniform')
        sca_bkg_choice = ones(20, 1)/20;
    elseif ~strcmp(sca_bkg_choice, 'rama') && ~isnumeric(sca_bkg_choice)
        error('script:scabkg', 'sca_bkg_choice should be ''rama'', ''uniform'', or a vector of 20 frequencies');
    end
end

% load the alignment
timer = tic;
alignment = loadfasta(['inputs/dhfr_' aln_choice '.fasta'], 'protein');

% load the PDB -- we'll use this repeatedly
pdb = pdbread('pdb/1RX2.pdb');
pdb_chain = 'A';

% get the Uniprot sequence to use as a reference for numbering
[up_seq, up_info] = getdbseq({pdb, pdb_chain});

% find the reference sequence in the alignment
[~, refidx] = seqsearch(alignment, 'annot', up_info.id, 'swaptofirst', false, ...
    'verbose', false);
if isempty(refidx)
    % if we can't find any matching annotation, don't give up: find the
    % best matching sequence (some alignments don't have annotations)
    alignment = seqsearch(alignment, 'seq', up_seq);
else
    % this swaps the reference sequence to the first position
    alignment = seqsearch(alignment, 'idx', refidx, 'seq', up_seq);
end
% truncate to the reference sequence and assign numbering
alignment = seqtruncate(alignment, 'seq', up_seq);
alignment.refseq(1).seqdb = 'uniprot';
alignment.refseq(1).seqid = up_info.accession;
disp(['Loading the alignment and truncating it took ' num2str(toc(timer)) ' seconds.']);

if filtering
    % we are also filtering the alignment
    timer = tic;
    % eliminate fragments
    alignment = filterrows(alignment, 'gaps', 0.2);
    % eliminate near-duplicates
    alignment = eliminatesimilarseqs(alignment, 0.8);
    % eliminate positions that most sequences don't have
    alignment = filtercolumns(alignment, 0.2);
    disp(['Filtering took ' num2str(toc(timer)) ' seconds.']);
end

% the sector from the original paper in Uniprot coordinates -- four versions,
% depending on the cutoff used
% (note that PDB coordinates are the same as Uniprot)
kimsecs_uniprot = {...
    [15 21 27 28 31 32 35 37 42 44 51 54 55 57 59 63 77 81 94 113 121 125 133], ...
    [11 15 21 22 27 28 31 32 35 37 39 40 42 44 47 50 51 54 55 56 57 59 63 64 77 81 94 100 111 113 121 122 125 126 133 153], ...
    [3 11 15 21 22 27 28 31 32 35 37 39 40 42 44 47 49 50 51 52 54 55 56 57 59 63 64 77 81 90 94 100 111 113 121 122 125 126 133 153], ...
    [3 5 11 13 15 18 21 22 25 27 28 31 32 35 37 39 40 42 44 45 47 49 50 51 52 53 54 55 56 57 59 63 64 77 81 90 94 95 99 100 107 111 113 121 122 125 126 133 153] ...
};

% get the PDB contact map
timer = tic;
coords = getpdbcoords(pdb, pdb_chain);
distmap = getdistmat(coords);
contactmap = (distmap < 6);

% restrict it to the same indices as the alignment
pdbmap = findseqmap(pdbgetseq(pdb, pdb_chain), up_seq, 'protein');
pdbidxs = findcommonidxs(pdbmap, alignment.refseq.map);
if length(pdbidxs) ~= size(alignment.data, 2)
    error('script:nopdb', 'Some alignment positions do not have PDB data.');
end
distmap = distmap(pdbidxs, pdbidxs);
contactmap = contactmap(pdbidxs, pdbidxs);
% we will use the atom coordinates, either restricted to the alignment or
% not, later
coords0 = coords;
coords = coords(pdbidxs);
disp(['Getting contact map from PDB took ' num2str(toc(timer)) ' seconds.']);

% a useful function
to_uniprot = @(secs) cellfun(@(sec) alignment.refseq.map(sec), secs, 'uniform', false);

% decide on a sector size as a number of residues
if sec_fraction < 1
    sec_count = round(size(alignment.data, 2)*sec_fraction);
else
    sec_count = sec_fraction;
end

%% SCA and background subtraction (if sub_rnd is true)

% SCA in Reynolds et al. uses the spectral norm and a thresholded version
% of the positional weights
if strcmp(sca_type, 'reynolds')
    % using SCA4.5, employed in Reynolds et al. 2011 -- reduction with spectral
    % norm, thresholded positional weights
    sca_options = {...
        'method', 'reduction', ...
        'redfct', @norm, ...
        'pwfct', @(f, patt) max(0, pwfctdkl(f, patt)) ...
        };
else
    sca_options = {'method', sca_type};
end
if ~strcmp(sca_bkg_choice, 'rama')
    % temporarily update the background distribution
    setdbbkg('protein', sca_bkg_choice);
end
sca0 = getsca(alignment, sca_options{:});

if sub_rnd
    timer = tic;
    % need default random number generator seed to get the exact same results
    rng('default');

    nrnd = 100;
    rndsca = cell(nrnd, 1);
    for i = 1:nrnd
        if mod(i, 10) == 0
            disp(['Working on sample ' int2str(i) ' of ' int2str(nrnd) '...']);
        end
        rndalign = alnscramble(alignment);
        rndsca{i} = getsca(rndalign, sca_options{:});
    end
    
    % subtract the mean random SCA matrix from the SCA matrix for the real
    % alignment
    Cavg = zeros(size(sca0.cmat));
    for i = 1:nrnd
        Cavg = Cavg + rndsca{i}.cmat;
    end
    Cavg = Cavg/nrnd;
    sca.cmat = sca0.cmat - Cavg;
    disp(['Subtraction of random SCA matrix expectation took ' num2str(toc(timer)) ' seconds.']);
else
    sca = sca0; %#ok<*UNRCH>
end

if ~strcmp(sca_bkg_choice, 'rama')
    % then revert to the usual one, to avoid confusions
    setdbbkg('protein', []);
end

%% get the sector

if ~use_top_ev
    % use the complicated method from Reynolds et al.
    kmax = 5;
    [ev, ~] = eigsorted(sca.cmat, kmax);
    
    % the cutoffs used in Reynolds et al. for the four sector sizes
    cutoffs = [0.0048 0.008 0.010 0.015];
    % XXX had to replace the first cutoff by 0.0048 (instead of 0.005) to
    % reproduce the Reynolds et al. residues exactly... this is probably
    % related to slightly different rounding and/or binning of the PDF in
    % getsectors vs. the Reynolds et al. code
    
    % Friedman-Diaconis rule
    binwidths = @(ev) 2*iqr(ev) * size(ev, 1)^(-0.33);
    % for each cutoff, find 'subsectors' corresponding to each of the top five
    % eigenvectors, then merge them together into a single sector
    % make this a function so we can reuse it later
    calculate_secs = @(ev, cutoffs) ...
        cellfun(...
        @(s) unique(cell2mat(s.secs)), ...
        arrayfun(...
        @(cutoff) getsectors(ev, ...
        'distrib', 'tlocationscale', ...
        'method', 'pdf', ...
        'tails', [1 2 2 2 2], ...
        'cutoff', cutoff ./ binwidths(ev)), ...
        cutoffs, 'uniform', false), ...
        'uniform', false);
    secs = calculate_secs(ev, cutoffs);
    % transform to Uniprot coordinates
    secs_uniprot = to_uniprot(secs);
    

    % this recovers the Reynolds et al. sector at all four cutoff levels exactly
    disp('Comparison between our sector and the one in Reynolds et al., at the four cutoff levels:');
    disp(arrayfun(@(i) setcompare(secs_uniprot{i}, kimsecs_uniprot{i}), 1:length(cutoffs)));
else
    % simply put a cutoff on the top eigenvector
    secopts = {...
        'method', 'count', ...
        'tails', 1, ...
        'cutoff', sec_count ...
    };
    secs0 = getsectors({sca, 1}, secopts{:});
    secs = secs0.secs;
    % transform to Uniprot coordinates
    secs_uniprot = to_uniprot(secs);

    % compare this to the Reynolds et al. sectors
    disp('Comparison between our sector and the one in Reynolds et al., at the four cutoff levels:');
    disp(arrayfun(@(i) setcompare(secs_uniprot{1}, kimsecs_uniprot{i}), 1:length(kimsecs_uniprot)));
end

%% 2. Get sectors of many sizes, for many cutoffs

if ~use_top_ev
    % this is a little more involved than one would expect because the sector
    % positions depend on the top five eigenvectors, so residues might overlap
    % between these, and thus it's not obvious how changing the cutoff affects
    % sector size
    
    % we use a simple-minded approach: we sweep the cutoff in some range,
    % searching for sectors of sizes that cover all the integers in the
    % corresponding interval
    timer = tic;
    % chose minimum and maximum cutoffs to use by trial-and-error
    cutoff1 = 1e-3;
    cutoff2 = 2e-2;
    smallsec = calculate_secs(ev, cutoff1);
    largesec = calculate_secs(ev, cutoff2);
    
    sizes = length(smallsec{1}):length(largesec{1});
    
    cutoffvec = zeros(size(sizes));
    cutoffvec(1) = cutoff1;
    cutoffvec(end) = cutoff2;
    
    manysecs = cell(size(sizes));
    manysecs{1} = smallsec{1};
    manysecs{end} = largesec{1};
    
    % keeping track of the number of times we called calculate_secs
    % minimum would be length(sizes) (if we knew which cutoffs to use)
    % we do a little over twice that, which isn't bad
    nevals = 2;
    for i = 1:length(sizes)
        crtsize = sizes(i);
        if ~isempty(manysecs{i})
            % we've calculated this already; go to the next one
            cutoff = cutoffvec(i);
            continue;
        end
        
        % binary search -- the sought-after sector size should correspond to a
        % cutoff between the one for the previous size and the maximum cutoff
        cutoff_left = cutoff;
        cutoff_right = cutoffvec(end);
        while true
            cutoff = (cutoff_left + cutoff_right) / 2;
            crtsec = calculate_secs(ev, cutoff);
            nevals = nevals + 1;
            size1 = length(crtsec{1});
            idx = size1 - sizes(1) + 1;
            if isempty(manysecs{idx})
                manysecs{idx} = crtsec{1};
                cutoffvec(idx) = cutoff;
            end
            if size1 < crtsize
                cutoff_left = cutoff;
            elseif size1 > crtsize
                cutoff_right = cutoff;
            else
                break;
            end
        end
    end
    
    %manysecs_uniprot = to_uniprot(manysecs);
    disp(['Calculation of sectors for many sizes took ' num2str(toc(timer)) '.']);
else
    sizes = 1:round(0.5*length(up_seq));
    manysecs = cell(size(sizes));
    for i = 1:length(sizes)
        crtsec = getsectors({sca, 1}, 'method', 'count', 'tails', 1, 'cutoff', sizes(i));
        manysecs{i} = crtsec.secs{1};
    end
end

%% Get conservations and conservation-based 'sectors' for many sizes

if ~strcmp(cons_bkg_choice, 'rama')
    % temporarily update the background distribution
    setdbbkg('protein', cons_bkg_choice);
end
cons = getconservation(alignment, 'method', 'kl', 'gaps', true);
if ~strcmp(cons_bkg_choice, 'rama')
    % temporarily update the background distribution
    setdbbkg('protein', []);
end

% helper function
findseclike = @(comps) arrayfun(@(s) ...
        getnth(1:s, getnth(2, @sort, comps, 'descend')), ...
    sizes, 'uniform', false);
% find conservation-based 'sectors'
manycons = findseclike(cons);

% find diag(SCA)-based 'sectors'
% make sure that we're taking the diagonal of the *unsubtracted* matrix
% otherwise diagsca is essentially zero
diagsca = diag(sca0.cmat);
manydiag = findseclike(diagsca);

%manycons_uniprot = to_uniprot(manycons);
%manydiag_uniprot = to_uniprot(manydiag);

% look, in particular, at the 'sectors' obtained from conservation and
% diag(SCA) that are the same size as the SCA sector
secscons = manycons(cellfun(@length, manysecs) == length(secs{1}));
secsdiag = manydiag(cellfun(@length, manysecs) == length(secs{1}));

%% Compare SCA sectors to those obtained from single-site statistics

simseccons = arrayfun(@(i) setcompare(manysecs{i}, manycons{i}), 1:length(manysecs));
simsecdiag = arrayfun(@(i) setcompare(manysecs{i}, manydiag{i}), 1:length(manysecs));

if make_figures
    figure;
    plot(sizes, simseccons, 'k', 'linewidth', 3);
    hold on;
    plot(sizes, simsecdiag, 'r', 'linewidth', 3);
    axis([0 max(sizes) + 4 0 1]);
    xlabel('Sector size');
    ylabel('Overlap with conservation or diag(SCA)-based ''sector''');
    legend(' Conservation', ' Diag(SCA)', 'location', 'best');
    legend(gca, 'boxoff');
    beautifygraph;
    set(gcf, 'name', 'Overlap between SCA sector and conservation or diag(SCA)-based one');
    preparegraph;
end

%% Check sector connectivity

[conn_secs, conninfo_secs] = checkconnectivity(secs{1}, contactmap, 'pvalue', true);
[conn_cons, conninfo_cons] = checkconnectivity(secscons{1}, contactmap, 'pvalue', true);
[conn_diag, conninfo_diag] = checkconnectivity(secsdiag{1}, contactmap, 'pvalue', true);

disp('Connectivity of SCA sector:');
disp([num2str(conn_secs, 2) ' (p = ' num2str(conninfo_secs.p, 2) ')']);
disp('Connectivity of conserved residues:');
disp([num2str(conn_cons, 2) ' (p = ' num2str(conninfo_cons.p, 2) ')']);
disp('Connectivity of residues with highest diag(SCA):');
disp([num2str(conn_diag, 2) ' (p = ' num2str(conninfo_diag.p, 2) ')']);

%% 3. Show dependence of conservation on position

residuepos = cell2mat(cellfun(@(s) mean(s.pos), coords, 'uniform', false));
cmpos = mean(residuepos);
resposcentered = bsxfun(@minus, residuepos, cmpos);
radii = arrayfun(@(i) norm(resposcentered(i, :)), 1:size(resposcentered, 1));

if make_figures
    figure;
    [~, ~, tmph] = scatterfit(radii, cons, ...
        'shapes', '.', 'colors', 'r', 'sizes', 20, ...
        'legend', 'c', 'legendloc', 'north', ...
        'corrtext', 'Correlation c = ', ...
        'showfit', false, 'showci', false);
    xlabel('distance from center');
    ylabel('conservation');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);
    preparegraph;
end

%% II. Comparing to experimental data

%% 1. Load the experimental data, and define some useful functions

% use a naming convention that only considers the first word in the header
expdata0 = csvload('inputs/dhfr_lightdep.csv', ...
    'namefct', @(s) lower(s(1:find(~isletter([s ' ']), 1)-1)));
expdata0_controls = csvload('inputs/dhfr_lightdep_controls.csv', ...
    'namefct', @(s) lower(s(1:find(~isletter([s ' ']), 1)-1)));

% select only the 'good' part of the data
mask = expdata0.data.good;
expdata = expdata0.data;
fnames = fieldnames(expdata0.data);
for i = 1:length(fnames)
    crtcol = expdata0.data.(fnames{i});
    expdata.(fnames{i}) = crtcol(mask);
end

mask = expdata0_controls.data.good;
fnames = fieldnames(expdata0_controls.data);
for i = 1:length(fnames)
    crtcol = expdata0_controls.data.(fnames{i});
    expdata_controls.(fnames{i}) = crtcol(mask);
end

% some results about the non-light-dependent controls
nld_mu = mean(expdata_controls.mean);
nld_std = std(expdata_controls.mean);
%nld_mu = 0.000648667;
%nld_std = 0.004658278;

% get a signed z-score
expdata.signedz = (expdata.mean - nld_mu) ./ sqrt((expdata.std.^2)/3 + nld_std^2);

% make a vector of the given size with ones at the given positions
makespikedvector0 = @(pos, n) full(sparse(pos(:), ones(size(pos(:))), ones(size(pos(:))), n, 1));
% make vector with size given by the size of the alignment
makespikedvector = @(pos) makespikedvector0(pos, size(alignment.data, 2));

%% 2. Get the list of 'touching' residues

% some helpful functions
ptinsph = @(p, c, r) norm(p - c) <= r;
manyptsinsph = @(res, c, r) arrayfun(@(i) ptinsph(res(i, :), c, r), 1:size(res, 1));
resinsph = @(res, c, r) any(manyptsinsph(res, c, r));
resinmanysph = @(res, cs, r) cellfun(@(c) resinsph(res, c, r), cs);
resinanysph = @(res, cs, r) any(resinmanysph(res, cs, r));
% taking advantage of the fact that PDB names are the same as the residue
% index in the PDB in this case
getcenters = @(pdbcoords, i, which) cellfun(...
    @(s) pdbcoords{i}.pos(strcmp(pdbcoords{i}.type, s), :), which, 'uniform', false);

timer = tic;
neighbors = cell(length(expdata.dhfr), 1);
for i = 1:length(expdata.dhfr)
    res2 = expdata.dhfr(i);
    res1 = res2 - 1;
    if res1 > 0
        centers = getcenters(coords0, res1, {'C', 'O'});
    else
        centers = {};
    end
    centers = [centers getcenters(coords0, res2, {'CA', 'N'})]; %#ok<AGROW>
    
    neighbors{i} = find(cellfun(@(res) resinanysph(res.pos, centers, cutoffrad), coords0));
end

expdata.residues = neighbors(:);
disp(['Finding the residue neighbors took ' num2str(toc(timer)) ' seconds.']);

%% 3. Recover the result from the paper

% first define functions for finding the residues 'touched' by a sector

% given a vector of weights for each position of an alignment, and a
% mapping between the alignment positions and PDB indices, calculate how
% much the "sector" identified by that vector of weights is "touching" the
% surface residues contained in the Reynolds et al. data

% more precisely: for each entry in the cell array 'residues' (each such
% entry should be a vector of integers corresponding to PDB indices of the
% positions touched by the corresponding surface residue), collect the
% values that correspond to those PDB positions from the vector v, and
% apply fct to them. return the resulting vector.
getsurfacetouching0 = @(v, ref, fct, residues) cellfun(...
    @(c) fct(v(ismember(ref, c))), residues);

% only the absolute value of the weights is used, and the weights are added
% when several alignment positions "touch" the same surface residue
getsurfacetouching = @(v) getsurfacetouching0(v, alignment.refseq.map, @(w) sum(abs(w)), expdata.residues);

% find all the surface positions that are touched by at least one sector
% residue
getfisher = @(touch_mask, sig_mask) [...
        sum(touch_mask & sig_mask) sum(~touch_mask & sig_mask) ; ...
        sum(touch_mask & ~sig_mask) sum(~touch_mask & ~sig_mask) ...
    ];
getfishertable = @(sec) getfisher(...
    (getsurfacetouching(makespikedvector(sec)) > 0), (expdata.z > zlimit));
getformatted = @(table) formattable(...
            [{'exp. sig.' ; 'not exp. sig.'} num2cell(table)], ...
            'header', {'', 'sector', 'non-sector'});
getstrings0 = @(ftable, table) {ftable, ...
    ['(Fisher p value: ' num2str(fishertest(table, 'tails', 1, 'type', 'abs'), 2) ...
            ', fraction sig. overall: ' num2str(sum(table(1, :)) / sum(table(:)), 2) ...
            ', fraction sig. in sector-touched residues: ' num2str(table(1, 1) / sum(table(:, 1)), 2) ')']};
getstrings1 = @(table) getstrings0(getformatted(table), table);
getstrings = @(sec) getstrings1(getfishertable(sec));

disp(' ');
disp('Contingency tables comparing SCA or conservation-based sectors to data');

disp('Sector:');
sectable = getfishertable(secs{1});
tmp = getstrings1(sectable);
disp(tmp{1});
disp(tmp{2});
disp(' ');

disp('Conservation:');
constable = getfishertable(secscons{1});
tmp = getstrings1(constable);
disp(tmp{1});
disp(tmp{2});
disp(' ');

disp('Diag(SCA):');
diagtable = getfishertable(secsdiag{1});
tmp = getstrings1(diagtable);
disp(tmp{1});
disp(tmp{2});
disp(' ');

disp('Chi^2 test comparing sector and conservation contingency tables: ');
disp(['p = ' num2str(chi2gof2([sectable(:) constable(:)]), 2)]);
disp(' ');

disp('Chi^2 test comparing sector and diag(SCA) contingency tables: ');
disp(['p = ' num2str(chi2gof2([sectable(:) diagtable(:)]), 2)]);
disp(' ');

%% 4. Calculate p-values for a range of sector sizes

secpvals = cellfun(@(sec) fishertest(getfishertable(sec), 'tails', 1, 'type', 'abs'), manysecs);
conspvals = cellfun(@(sec) fishertest(getfishertable(sec), 'tails', 1, 'type', 'abs'), manycons);
diagpvals = cellfun(@(sec) fishertest(getfishertable(sec), 'tails', 1, 'type', 'abs'), manydiag);

getchi2 = @(table1, table2) chi2gof2([table1(:) table2(:)]);
[comppvals, compchi2] = arrayfun(@(i) getchi2(getfishertable(manysecs{i}), getfishertable(manycons{i})), 1:length(manysecs));
[diagcomppvals, diagcompchi2] = arrayfun(@(i) getchi2(getfishertable(manysecs{i}), getfishertable(manydiag{i})), 1:length(manysecs));

%% make figures for chi^2 p-values for sectors vs. conservation

if make_figures
	figure;
    plot(sizes, comppvals, '.r', 'markersize', 20);
    tmp = axis;
    axis([0 max(sizes) + 4 0 floor(tmp(4)*10)/10]);
    xlabel('Sector size');
    ylabel('p-value');
    set(gca, 'xgrid', 'on', 'ygrid', 'on');
    beautifygraph;
    set(gcf, 'name', 'Comparison of sector to conservation using a chi^2 test');
    preparegraph;
end

%% make figures for chi^2 p-values for sectors vs. diag(SCA)

if make_figures
	figure;
    plot(sizes, diagcomppvals, '.r', 'markersize', 20);
    tmp = axis;
    axis([0 max(sizes) + 4 0 floor(tmp(4)*10)/10]);
    xlabel('Sector size');
    ylabel('p-value');
    set(gca, 'xgrid', 'on', 'ygrid', 'on');
    beautifygraph;
    set(gcf, 'name', 'Comparison of sector to diag(SCA) using a chi^2 test');
    preparegraph;
end