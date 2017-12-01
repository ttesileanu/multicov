% This script performs the analysis of the SCA method on the PDZ
% alignments, as described in the paper. This is also the script used to
% generate the graphs.
%
% The basic purpose of this script is to show how the sector found in
% McLaughlin Jr. et al. 2012 compares to the top conserved positions, how
% the comparison changes if we change the alignment or the threshold, etc.
%
% There are some options that can be set by defining some variables before
% running section I.0 of the script. Leave these variables undefined to get
% the default behavior.
%   'aln_choice'
%       Choice of alignment to use. Can be
%         'rama': the original alignment from McLaughlin et al. 2012,
%         'hhblits_e5': use HHblits alignment with e5 threshold,
%         'hhblits_e10': use HHblits alignment with e10 threshold,
%         'pfam': raw Pfam alignment.
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
%   'data_choice'
%       Which ligand to use.
%         'cript': the cognate ligand to PSD95
%         'tm2f':  the mutated ligand
%         'diff':  diference between the datasets for  the two ligands
%       (default: 'cript')
%   'data_restrict'
%       Restrict the data to only a particular kind of mutations. This is a
%       character string describing the kinds of mutants to look at (that
%       is, the allowed target amino acids).
%       (default: no restriction)
%   'filtering'
%       Whether to filter the alignment. This implies first capping the
%       fraction of gaps per sequence to 20%, then eliminating sequences
%       that are more than 80% identical, and finally removing columns of
%       the alignment that contain more than 20% gaps.
%       (default: false)
%   'make_figures'
%       Whether to make the figures.
%       (default: true)
%   'nfctpos'
%       Number of positions to consider functionally significant.
%       (default: 20)
%   'sca_bkg_choice'
%       A choice of background distribution to be used with SCA. See
%       'cons_bkg_choice' for choices.
%       (default: 'rama')
%   'sca_type'
%       Type of SCA matrix to use. Options:
%        'reduction':   reduction method SCA. (version 4.0)
%        'projection':  projection method SCA. (version 5.0)
%        'binary':      binary approximation SCA.
%       (default: 'projection')
%   'sec_fraction'
%       Fraction of protein to be used for the sector. The protein size is
%       assumed to be equal to the width of the alignment. If sec_fraction
%       is set to a number larger than or equal to 1, it is assumed to be
%       equal to the number of residues that should be taken up by the
%       sector.
%       (default: 0.25)
%   'use_rama_ats'
%       If aln_choice is 'rama' and this option is set to true, we use a
%       slightly different filtering convention to match the exact
%       alignment used in the paper.
%       (default: true)

% Tiberiu Tesileanu (2013-2014)

%% I. Information we get from the alignment alone

%% 0. Setup

% handle some options
if ~exist('aln_choice', 'var')
    aln_choice = 'hhblits_e10';
end
if ~exist('data_choice', 'var')
    data_choice = 'cript';
end
if ~exist('data_restrict', 'var')
    data_restrict = alphagetletters('protein', 'nogap');
end
if ~exist('filtering', 'var')
    filtering = false;
end
if ~exist('use_rama_ats', 'var')
    use_rama_ats = true;
end
if ~exist('make_figures', 'var')
    make_figures = true;
end
if ~exist('nfctpos', 'var')
    nfctpos = 20;
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
if ~exist('sca_type', 'var')
    sca_type = 'projection';
end
if ~exist('sec_fraction', 'var')
    sec_fraction = 0.25;
end

% load the alignment
timer = tic;
alignment = loadfasta(['inputs/pdz_' aln_choice '.fasta'], 'protein');
%alignment = loadfasta(['inputs/pdz_' aln_choice '.fasta'], 'protein', ...
%    'keepmask', @(v) true(size(v)), 'postprocess', @(v) upper(v));

% we might need a copy of the original alignment
alignment0 = alignment;

% load the PDB -- we'll use this repeatedly
pdb = pdbread('pdb/1BE9.pdb');
pdb_chain = 'A';

% get the Uniprot sequence to use as a reference for numbering
[up_seq, up_info] = getdbseq({pdb, pdb_chain});
% but! make sure to search only the subset of the sequence that is relevant
% for the PDB structure (some Uniprot sequences are much longer and can
% confuse the search -- for PDZ there are 3 different sequences in the Pfam
% alignment that match different portions of the Uniprot DLG4_RAT sequence;
% only one of these is contained in the PDB)
up_seq_trunc = up_seq(up_info.range);

% find the reference sequence in the alignment
[~, refidx] = seqsearch(alignment, 'annot', up_info.id, 'swaptofirst', false, ...
    'verbose', false);
if isempty(refidx)
    % if we can't find any matching annotation, don't give up: find the
    % best matching sequence (some alignments don't have annotations)
    alignment = seqsearch(alignment, 'seq', up_seq_trunc);
else
    % this swaps the reference sequence to the first position
    alignment = seqsearch(alignment, 'idx', refidx, 'seq', up_seq_trunc);
end
% truncate to the reference sequence and assign numbering
alignment = seqtruncate(alignment, 'seq', up_seq_trunc, 'posnames', up_info.range);
% update database information
alignment.refseq(1).seqdb = 'uniprot';
alignment.refseq(1).seqid = up_info.accession;
disp(['Loading alignment and truncating to PDB positions took ' num2str(toc(timer)) ' seconds.']);

if strcmp(aln_choice, 'rama') && use_rama_ats
    % in McLaughlin et al. they use a more complex method for matching
    % alignment columns to PDB... the result of that method is that one
    % column gets a mismatched PDB coordinate. We reproduce that here, for
    % consistency with the published work. The results are not significantly
    % affected by this choice.
    
    % first, McLaughlin et al. eliminate columns that have too many gaps
    alignment = filtercolumns(alignment, 0.2);
    
    % manually add in the misidentified position
    newdata = repmat(blanks(size(alignment.data, 2) + 1), size(alignment.data, 1), 1);
    newdata(:, 1:30) = alignment.data(:, 1:30);
    newdata(:, 32:end) = alignment.data(:, 31:end);
    newdata(:, 31) = alignment0.data(:, 51);
    % alignment0 had the focal sequence in a different position; seqsearch
    % swapped it in first position. we now need to do this swap manually
    % for column 31
    refidx = 77;
    newdata([1 refidx], 31) = newdata([refidx 1], 31);
    
    alignment.alphawidths = size(newdata, 2);
    alignment.data = newdata;
    alignment.refseq.map = [alignment.refseq.map(1:30) ; 334 ; alignment.refseq.map(31:end)];
elseif filtering
    % here we are also filtering the alignment
    timer = tic; %#ok<*UNRCH>
    % eliminate fragments
    alignment = filterrows(alignment, 'gaps', 0.2);
    % eliminate near-duplicates
    alignment = eliminatesimilarseqs(alignment, 0.8);
    % eliminate positions that most sequences don't have
    alignment = filtercolumns(alignment, 0.2);
    disp(['Filtering took ' num2str(toc(timer)) ' seconds.']);
end

% get the PDB contact map
timer = tic;
coords = getpdbcoords(pdb, pdb_chain);
distmap = getdistmat(coords);
contactmap = (distmap < 6);

% restrict it to the same indices as the alignment
pdbmap_trunc = findseqmap(pdbgetseq(pdb, pdb_chain), up_seq_trunc, 'protein');
pdbmap = zeros(size(pdbmap_trunc));
pdbmap(pdbmap_trunc ~= 0) = up_info.range(pdbmap_trunc(pdbmap_trunc ~= 0));
pdbidxs = findcommonidxs(pdbmap, alignment.refseq.map);
if length(pdbidxs) ~= size(alignment.data, 2)
    error('script:nopdb', 'Some alignment positions do not have PDB data.');
end
distmap = distmap(pdbidxs, pdbidxs);
contactmap = contactmap(pdbidxs, pdbidxs);
disp(['Loading PDB file and getting contact map took ' num2str(toc(timer)) ' seconds.']);

% decide on a sector size as a number of residues
if sec_fraction < 1
    sec_count = round(size(alignment.data, 2)*sec_fraction);
else
    sec_count = sec_fraction;
end

%% 1. Reproducing the sector from McLaughlin Jr. et al.

% the sector from the original paper in Uniprot/PDB coordinates
ramasec_uniprot = [322 ; 323 ; 325 ; 327 ; 329 ; 330 ; 336 ; 347 ; 351 ; 353 ; 359 ; 362 ; 363 ; 364 ; 372 ; 375 ; 376 ; 379 ; 386 ; 388];

% run SCA on the alignment
if ~strcmp(sca_bkg_choice, 'rama')
    % temporarily update the background distribution
    setdbbkg('protein', sca_bkg_choice);
end
sca = getsca(alignment, 'method', sca_type);
if ~strcmp(sca_bkg_choice, 'rama')
    % then revert to the usual one, to avoid confusions
    setdbbkg('protein', []);
end
[ev, ~] = eigsorted(sca.cmat);

% we will use this to calculate sector-like objects not only from SCA, but
% also from single-site statistics
calculate_secs = @(ev) ...
    getsectors(ev(:, 1), ...
    'method', 'count', ...
    'distrib', 'none', ...
    'cutoff', sec_count, ...
    'tails', 1, ...
    'refseq', alignment.refseq);

secssca = calculate_secs(ev);

% this reproduces the sector from the paper exactly when aln_choice is 'rama'
% and use_rama_ats is true
disp(['Identity between sector obtained here and sector from McLaughlin Jr. et al.: ' num2str(comparesectors(secssca.secs_mapped{1}, ramasec_uniprot)) '.']);

%% 2. Relation between top SCA eigenvector and conservation

timer = tic;
if ~strcmp(cons_bkg_choice, 'rama')
    % temporarily update the background distribution
    setdbbkg('protein', cons_bkg_choice);
end
cons = getconservation(alignment, 'method', 'kl', 'gaps', true);
cons_max = getconservation(alignment, 'method', 'max');
if ~strcmp(cons_bkg_choice, 'rama')
    % temporarily update the background distribution
    setdbbkg('protein', []);
end
% find a conservation-based 'sector'
secscons = calculate_secs(cons);
% find a diag(SCA)-based 'sector'
secsdiag = calculate_secs(diag(sca.cmat));
disp(['Conservation calculations took ' num2str(toc(timer)) ' seconds.']);

% compare sectors obtained from SCA to conservation / diag(SCA)
disp('Comparison of SCA sector with conservation-based sector:');
disp(comparesectors(secssca.secs, secscons.secs));

disp('Comparison of SCA sector with diag(SCA)-based sector:');
disp(comparesectors(secssca.secs, secsdiag.secs));

%% compare different definitions for conservation

if make_figures
    figure;
    [~, ~, tmph] = scatterfit(cons, cons_max, ...
        'shapes', '.', ...
        'colors', 'r', ...
        'sizes', 20, ...
        'legend', 'c', ...
        'corrtext', 'Correlation c = ', ...
        'showci', false, ...
        'style', {'k--', 'linewidth', 2}, ...
        'legendloc', 'north');
    xlabel('Conservation (KL divergence)');
    ylabel('Conservation (frequency of consensus amino acid)');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);
    set(gcf, 'name', 'Different definitions for conservation');
    preparegraph;
end

%% compare top eigenvector of SCA to single-site statistics

if make_figures
    % show the correlation between the top eigenvector and conservation
    figure;
    [~, ~, tmph] = ...
        scatterfit(cons, ev(:, 1), ...
        'shapes', '.', ...
        'colors', 'r', ...
        'sizes', 20, ...
        'legend', 'c', ...
        'corrtext', 'Correlation c = ', ...
        'showfit', false, 'showci', false, ...
        'legendloc', 'north');
    tmp = axis;
    axis([0 tmp(2) 0 tmp(4)]);
    xlabel('Conservation');
    ylabel('Components of top eigenvector');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);
    set(gcf, 'name', 'Conservation vs. top eigenvector of SCA matrix');
    preparegraph;
    
    % show the correlation between the top eigenvector and diag(SCA)
    figure;
    [~, ~, tmph] = ...
        scatterfit(sqrt(diag(sca.cmat)), ev(:, 1), ...
        'shapes', '.', ...
        'colors', 'r', ...
        'sizes', 20, ...
        'legend', 'c', ...
        'corrtext', 'Correlation c = ', ...
        'showfit', false, 'showci', false, ...
        'legendloc', 'north');
    tmp = axis;
    axis([0 tmp(2) 0 tmp(4)]);
    xlabel('Square root of diagonal of SCA matrix');
    ylabel('Components of top eigenvector');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);
    set(gcf, 'name', 'Diagonal vs. top eigenvector of SCA matrix');
    preparegraph;
    
    % show the correlation between diag(SCA) and conservation
    figure;
    [~, ~, tmph] = ...
        scatterfit(cons, sqrt(diag(sca.cmat)), ...
        'shapes', '.', ...
        'colors', 'r', ...
        'sizes', 20, ...
        'legend', 'c', ...
        'corrtext', 'Correlation c = ', ...
        'showfit', false, 'showci', false, ...
        'legendloc', 'north');
    tmp = axis;
    axis([0 tmp(2) 0 tmp(4)]);
    xlabel('Conservation');
    ylabel('Square root of diagonal elements');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);
    set(gcf, 'name', 'Conservation vs. diagonal of SCA matrix');
    preparegraph;
end

%% 3. Check the structural connectivity of the sector-like objects

[conn_secs, conninfo_secs] = checkconnectivity(secssca.secs{1}, contactmap, 'pvalue', true);
[conn_cons, conninfo_cons] = checkconnectivity(secscons.secs{1}, contactmap, 'pvalue', true);
[conn_diag, conninfo_diag] = checkconnectivity(secsdiag.secs{1}, contactmap, 'pvalue', true);

disp('Connectivity of SCA sector:');
disp([num2str(conn_secs, 2) ' (p = ' num2str(conninfo_secs.p, 2) ')']);
disp('Connectivity of conserved residues:');
disp([num2str(conn_cons, 2) ' (p = ' num2str(conninfo_cons.p, 2) ')']);
disp('Connectivity of residues with highest diag(SCA):');
disp([num2str(conn_diag, 2) ' (p = ' num2str(conninfo_diag.p, 2) ')']);

%% II. Comparing to experimental binding affinity data

%% 1. Load the data

expresults = open('inputs/McLaughlin_Ranganathan_data.mat');

% calculate per-site averages
aamask_cript = ismember(expresults.single_mutagenesis_cript.aa, data_restrict);
aamask_Tm2F = ismember(expresults.single_mutagenesis_Tm2F.aa, data_restrict);
expresults.single_mutagenesis_cript.Eavgx = mean(expresults.single_mutagenesis_cript.data(aamask_cript, :), 1);
expresults.single_mutagenesis_Tm2F.Eavgx = mean(expresults.single_mutagenesis_Tm2F.data(aamask_cript, :), 1);

% select the data source
switch data_choice
    case 'cript'
        expdata = expresults.single_mutagenesis_cript;
    case 'tm2f'
        expdata = expresults.single_mutagenesis_Tm2F;
    case 'diff'
        expdata = expresults.single_mutagenesis_Tm2F;
        expdata.data = expdata.data - expresults.single_mutagenesis_cript;
        expdata.Eavgx = expdata.Eavgx - expresults.single_mutagenesis_cript.Eavgx;
    otherwise
        error(['Unrecognized data choice ' data_choice '.']);
end

%% 2. Compare experimental data to SCA top eigenvector and conservation

if make_figures
    figure;
    
    set(gcf, 'position', [100 200 800 300]);
    subplot(121);
    scatterfit(ev(:, 1), expdata.Eavgx, ...
        'shapes', '.', 'colors', 'r', 'sizes', 20, ...
        'style', {'k--', 'linewidth', 2}, ...
        'refseq', {alignment.refseq.map, expdata.pos}, ...
        'legend', 'csp');
    xlabel('Component on top SCA eigenvector');
    ylabel(['Experimental mutational effect (' data_choice ')']);
    beautifygraph;
    
    subplot(122);
    scatterfit(cons, expdata.Eavgx, ...
        'shapes', '.', 'colors', 'r', 'sizes', 20, ...
        'style', {'k--', 'linewidth', 2}, ...
        'refseq', {alignment.refseq.map, expdata.pos}, ...
        'legend', 'csp');
    xlabel('Conservation');
    ylabel(['Experimental mutational effect (' data_choice ')']);
    beautifygraph;
    set(gcf, 'name', 'Experimental data vs. conservation and top eigenvector SCA');
    preparegraph;
end

%% 3. Compare SCA and conservation to data using Fisher's exact test

% set up some useful functions

% find set of "functionally significant" positions -- e.g., top 20 ones, as
% in the paper
fctidxs = getnth(1:nfctpos, getnth(2, @sort, expdata.Eavgx));
fctpdb0 = intersect(expdata.pos(fctidxs), alignment.refseq.map);
% total number of positions that appear both in the alignment and in the
% experimental data
ntotpos0 = length(intersect(expdata.pos, alignment.refseq.map));

% keep only residues that also have experimental data
keepresidueswithexp = @(sec) intersect(sec, expdata.pos);
% helper function
tablefromthree = @(a, b, c, n) [a, b ; c, n - a - b - c];
generatetable0 = @(secwithexp) tablefromthree(...
    length(intersect(secwithexp, fctpdb0)), length(setdiff(fctpdb0, secwithexp)), ...
    length(setdiff(secwithexp, fctpdb0)), ntotpos0);
% generate a contingency table for comparing sector to functional residues
generatetable = @(sec) generatetable0(keepresidueswithexp(sec));
% generate message showing Fisher p value
generatepvaltext = @(ftable) ['(Fisher p value: ' num2str(fishertest(ftable, 'tails', 1, 'type', 'abs'), 2) ...
            ', fraction sig. overall: ' num2str(sum(ftable(1, :)) / sum(ftable(:)), 2) ...
            ', fraction sig. in sector residues: ' num2str(ftable(1, 1) / sum(ftable(:, 1)), 2) ')'];
% generate formatted table
generateformattedtab = @(ftable) formattable(...
    [{'exp. sig.' ; 'not exp. sig.'} num2cell(ftable)], ...
    'header', {'', 'sector', 'non-sector'});

disp(' ');
disp('Contingency tables comparing SCA or conservation-based sectors to data:');

% see how well the sector does
disp('Sector:');
secftable = generatetable(secssca.secs_mapped{1}{1});
disp(generateformattedtab(secftable));
disp(generatepvaltext(secftable));
disp(' ');

% make a Fisher's exact table/test for conservation
disp('Conservation:');
consftable = generatetable(secscons.secs_mapped{1}{1});
disp(generateformattedtab(consftable));
disp(generatepvaltext(consftable));
disp(' ');
        
% make a Fisher's exact table/test for diag(SCA)
disp('Diagonal of SCA matrix:');
diagftable = generatetable(secsdiag.secs_mapped{1}{1});
disp(generateformattedtab(diagftable));
disp(generatepvaltext(diagftable));

% Compare sector to conserved residues
disp(' ');
disp('Chi^2 p-value for comparing sector to conservation contingency tables:');
disp(['p = ' num2str(chi2gof2([secftable(:) consftable(:)]), 2)]);

disp(' ');
disp('Chi^2 p-value for comparing sector to diag(SCA) contingency tables:');
disp(['p = ' num2str(chi2gof2([secftable(:) diagftable(:)]), 2)]);

%% 4. Histogram of mutational effects for sector vs. all residues

% recreating (parts of) Figure 3 of McLaughlin et al.
disp(' ');
disp('Histogram comparisons between SCA and conservation-based sectors and data');
pdzshowhistcomp(expdata, alignment.refseq.map, ev(:, 1), cons, [], 'pvals', true, ...
    'format', '%.1e', 'figures', make_figures, ...
    'n', length(secssca.secs{1}));

% what happens if we change the number of residues?
% if make_figures
%     pdzshowhistcomp(expdata, alignment.refseq.map, evproj(:, 1), cons, [], 'n', 15, 'pvals', true, 'format', '%.1e');
%     pdzshowhistcomp(expdata, alignment.refseq.map, evproj(:, 1), cons, [], 'n', 25, 'pvals', true, 'format', '%.1e');
%     pdzshowhistcomp(expdata, alignment.refseq.map, evproj(:, 1), cons, [], 'n', 30, 'pvals', true, 'format', '%.1e');
% end

%% 5. Perform the comparison for a whole range of sector sizes

% range of sizes to use
nmin = 1;
nmax = 50;

[idxs0_aln, idxs0_exp] = findcommonidxs(alignment.refseq.map, expdata.pos);
values0_exp = expdata.Eavgx(idxs0_exp);

count = nmax - nmin + 1;
sec_ks = zeros(1, count);
sec_wwp = zeros(1, count);
cons_ks = zeros(1, count);
cons_wwp = zeros(1, count);
comparison_ks = zeros(1, count);
comparison_wwp = zeros(1, count);
sec_pval = zeros(1, count);
cons_pval = zeros(1, count);
chi2_pval = zeros(1, count);

sec_wwp_stats = cell(1, count);
cons_wwp_stats = cell(1, count);
comparison_wwp_stats = cell(1, count);
chi2_stat = zeros(1, count);
for n = nmin:nmax
    [~, sorted_vec] = sort(ev(:, 1), 'descend');
    selected_sec = sorted_vec(1:n);
    values0_sec = values0_exp(ismember(idxs0_aln, selected_sec));
    
    [~, sorted_vec] = sort(cons(:, 1), 'descend');
    selected_cons = sorted_vec(1:n);
    values0_cons = values0_exp(ismember(idxs0_aln, selected_cons));
    
    i = n - nmin + 1;
    [~, comparison_ks(i)] = kstest2(values0_sec, values0_cons, 0.05);
    [comparison_wwp(i), ~, comparison_wwp_stats{i}] = ranksum(values0_sec, values0_cons, 'method', 'approximate');
    
    [~, sec_ks(i)] = kstest2(values0_sec, values0_exp, 0.05, 'larger');
    [~, cons_ks(i)] = kstest2(values0_cons, values0_exp, 0.05, 'larger');
    
    [sec_wwp(i), ~, sec_wwp_stats{i}] = ranksum(values0_sec, values0_exp, 'tail', 'left');
    [cons_wwp(i), ~, cons_wwp_stats{i}] = ranksum(values0_cons, values0_exp, 'tail', 'left');
    
    crtsec_mapped = alignment.refseq.map(selected_sec);
    crtcons_mapped = alignment.refseq.map(selected_cons);
    
    secftable = generatetable(crtsec_mapped);
    consftable = generatetable(crtcons_mapped);
    sec_pval(i) = fishertest(secftable, 'tails', 1, 'type', 'abs');
    cons_pval(i) = fishertest(consftable, 'tails', 1, 'type', 'abs');
    
    [chi2_pval(i), chi2_stat(i)] = chi2gof2([secftable(:) consftable(:)]);
end

%% make plot of Mann Whitney p-values

if make_figures
    figure;
    plot(nmin:nmax, comparison_wwp, '.r', 'markersize', 20);
    axis([nmin nmax 0 1]);
    xlabel('Sector size');
    ylabel('p-value');
    legend(gca, 'boxoff');
    set(gca, 'xgrid', 'on', 'ygrid', 'on');
    beautifygraph;
    set(gcf, 'name', 'Sector vs. conserved comparison (Mann Whitney U test)');
    preparegraph;
end

%% make plot of ranksum statistic

if make_figures
    figure;
    plot(nmin:nmax, cellfun(@(s) s.zval, comparison_wwp_stats), '.k', 'markersize', 20);
    hold on;
    xlabel('Sector size');
    ylabel('Rank sum');
    tmph = fill([nmin+0.2 nmax-0.2 nmax-0.2 nmin+0.2], [2 2 -2 -2], [1 0.7 0.7], 'edgecolor', 'none');
    tmph1 = plot([nmin nmax], [2 2], 'r', 'linewidth', 2);
    tmph2 = plot([nmin nmax], [-2 -2], 'r', 'linewidth', 2);
    tmph3 = plot([nmin nmax], [0 0], 'k');
    set(gca, 'xgrid', 'on', 'ygrid', 'on');
    uistack(tmph3, 'bottom');
    uistack([tmph1 tmph2], 'bottom');
    uistack(tmph, 'bottom');
    axis([nmin nmax -4 4]);
    beautifygraph;
    set(gcf, 'name', 'Sector vs. conserved comparison (Mann Whitney U test)');
    preparegraph;
end

%% 6. Look at other eigenvectors

[idxs0_aln, idxs0_exp] = findcommonidxs(alignment.refseq.map, expdata.pos);

% get those portions of the eigenvectors that match the available
% mutagenesis data
ev_masked = ev(idxs0_aln, :);
cript_masked = expresults.single_mutagenesis_cript.Eavgx(idxs0_exp);
tm2f_masked = expresults.single_mutagenesis_Tm2F.Eavgx(idxs0_exp);

% we use the top three eigenvectors as predictors of the two binding
% affinities (to CRIPT and to Tm2F)
predictors = [ev_masked(:, 1:3) ones(size(ev_masked, 1), 1)];
[fit_cript, int_cript, ~, ~, stats_cript] = regress(cript_masked(:), predictors);
[fit_tm2f, int_tm2f, ~, ~, stats_tm2f] = regress(tm2f_masked(:), predictors);

% some helper functions for displaying the plot labels below
signchars = '-+';
choosesgn = @(x) signchars(1 + ceil((sign(x) + 1)/2));
stars = {'*', ''};
starchoose = @(name) stars{isempty(name) + 1};
lincombtostr = @(fdata, names) [num2str(fdata(1), 3) '*' names{1} ' ' ...
    cell2mat(arrayfun(@(i) [choosesgn(fdata(i)) ' ' num2str(abs(fdata(i)), 3) starchoose(names{i}) names{i} ' '], 2:length(fdata), 'uniform', false))];

%% plots

if make_figures
    figure;
    [~, ~, tmph] = scatterfit(...
        predictors*fit_cript, cript_masked, ...
        'shapes', '.', 'colors', 'r', ...
        'sizes', 20, 'legend', 'c', ...
        'corrtext', 'Correlation c = ', ...
        'legendloc', 'north', ...
        'showfit', true, 'showci', false, ...
        'line', [1 0], 'style', {'--k'});
    xlabel(lincombtostr(fit_cript, {'v1', 'v2', 'v3', ''}));
    ylabel('Experimental effect (CRIPT)');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);
    tmp = get(gcf, 'position');
    set(gcf, 'position', [tmp(1:2) tmp(3:4)*0.75]);
    preparegraph;
    
    figure;
    [~, ~, tmph] = scatterfit(...
        predictors*fit_tm2f, tm2f_masked, ...
        'shapes', '.', 'colors', 'r', ...
        'sizes', 20, 'legend', 'c', ...
        'corrtext', 'Correlation c = ', ...
        'legendloc', 'north', ...
        'showfit', true, 'showci', false, ...
        'line', [1 0], 'style', {'--k'});
    xlabel(lincombtostr(fit_tm2f, {'v1', 'v2', 'v3', ''}));
    ylabel('Experimental effect (Tm2F)');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);
    tmp = get(gcf, 'position');
    set(gcf, 'position', [tmp(1:2) tmp(3:4)*0.75]);
    preparegraph;
end