% This script performs the analysis of the SCA method on a LacI
% alignment, as described in the paper. This is also the script used to
% generate the graphs.
%
% The basic purpose of this script is to show how the sector positions
% and the top conserved positions compare to the mutagenesis data from
% Markiewicz et al. 1994.
%
% There are some options that can be set by defining some variables before
% running the first section of the script. Leave these variables undefined
% to get the default behavior.
%   'aln_choice'
%       Choice of alignment to use. Can be
%         'hhblits_e10': use HHblits alignment with e10 threshold,
%         'hhblits_e30': use HHblits alignment with e30 threshold.
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
%   'data_restrict'
%       Restrict the data to only a particular kind of mutation. This is a
%       character string describing the kinds of mutants to look at (that
%       is, the allowed target amino acids). For example, set data_restrict
%       to 'A' to focus on the alanine scan part of the data.
%       (default: no restriction)
%   'fctthreshold'
%       How many amino acid substitutions need to be deleterious for a site
%       to be considered functionally-significant.
%       (default: min(8, length(data_restrict)))
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
%       (default: 'projection')
%   'sec_fraction'
%       Fraction of protein to be used for the sector. The protein size is
%       assumed to be equal to the width of the alignment. If sec_fraction
%       is set to a number larger than or equal to 1, it is assumed to be
%       equal to the number of residues that should be taken up by the
%       sector.
%       (default: 0.25)

% Tiberiu Tesileanu (2014)

%% Setup

if ~exist('aln_choice', 'var')
    aln_choice = 'hhblits_e10';
end
if ~exist('filtering', 'var')
    filtering = false;
end
if ~exist('sca_type', 'var')
    sca_type = 'projection';
end
if ~exist('sec_fraction', 'var')
    sec_fraction = 0.25;
end
if ~exist('make_figures', 'var')
    make_figures = true;
end
if ~exist('data_restrict', 'var')
    data_restrict = alphagetletters('protein', 'nogap');
end
if ~exist('fctthreshold', 'var')
    fctthreshold = min(8, length(data_restrict));
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

timer = tic;
alignment = loadfasta(['inputs/laci_' aln_choice '.fasta'], 'protein');
%alignment = loadfasta(['inputs/laci_' aln_choice '.fasta'], 'protein', ...
%    'keepmask', @(v) true(size(v)), 'postprocess', @(v) upper(v));
disp(['Loading the alignment took ' num2str(toc(timer)) ' seconds.']);

timer = tic;
pdb = pdbread('pdb/1JWL.pdb');
pdb_chain = 'A';
% get the Uniprot sequence to use as a reference for numbering
[up_seq, up_info] = getdbseq({pdb, pdb_chain});

% keep only the residues for which there are mutants in Markiewicz et al.
datarange = 2:329;
up_seq_trunc = up_seq(datarange);

% place LACI_ECOLI sequence in first position and assign mapping
alignment = seqsearch(alignment, 'annot', up_info.id, 'seq', up_seq_trunc);
alignment = seqtruncate(alignment, 'seq', up_seq_trunc, 'posnames', datarange);
% update database information
alignment.refseq(1).seqdb = 'uniprot';
alignment.refseq(1).seqid = up_info.accession;

% make a truncated alignment if required
if filtering
    % this takes a rather long time, so if this is run multiple times, it's
    % a good idea to do it once, save it, and then load from file
%    load('generated/laci_hhblits_e30_filtered.mat'); %#ok<*UNRCH>
    timer2 = tic; %#ok<*UNRCH>
    % eliminate fragments
    alignment = filterrows(alignment, 'gaps', 0.2);
    % eliminate near-duplicates
    alignment = eliminatesimilarseqs(alignment, 0.8);
    % eliminate positions that most sequences don't have
    alignment = filtercolumns(alignment, 0.2);

%    save('generated/laci_hhblits_e30_filtered.mat', 'alignment');
    disp(['Filtering took ' num2str(toc(timer2)) ' seconds.']);
end

% get the PDB contact map
coords = getpdbcoords(pdb, pdb_chain);
distmap = getdistmat(coords);
contactmap = (distmap < 6);

% restrict it to the same indices as the alignment
pdbmap = findseqmap(pdbgetseq(pdb, pdb_chain), up_seq, 'protein');
pdbidxs = findcommonidxs(pdbmap, alignment.refseq.map);
if length(pdbidxs) ~= size(alignment.data, 2)
    warning('script:nopdb', 'Some alignment positions do not have PDB data. Connectivity information may be wrong.');
end
distmap = distmap(pdbidxs, pdbidxs);
contactmap = contactmap(pdbidxs, pdbidxs);

binalign = alntobin(alignment);
disp(['Truncating, filtering alignment, and getting contact map took ' num2str(toc(timer)) ' seconds.']);

% decide on a sector size as a number of residues
if sec_fraction < 1
    sec_count = round(size(alignment.data, 2)*sec_fraction);
else
    sec_count = sec_fraction;
end

%% What does SCA make of this?

timer = tic;
if ~strcmp(sca_bkg_choice, 'rama')
    % temporarily update the background distribution
    setdbbkg('protein', sca_bkg_choice);
end
sca = getsca(binalign, 'method', sca_type);
if ~strcmp(sca_bkg_choice, 'rama')
    % then revert to the usual one, to avoid confusions
    setdbbkg('protein', []);
end
[ev, lb] = eigsorted(sca.cmat);

% Get the sectors
secopts = {...
    'method', 'count', ...
    'tails', 1, ...
    'cutoff', sec_count ...
};
secs = getsectors({sca, 1}, secopts{:});
disp(['Calculating the SCA matrix, its spectrum, and the sector took ' num2str(toc(timer)) ' seconds.']);

%% Get conservation

timer = tic;
if ~strcmp(cons_bkg_choice, 'rama')
    % temporarily update the background distribution
    setdbbkg('protein', cons_bkg_choice);
end
cons = getconservation(binalign, 'method', 'kl', 'gaps', true);
cons_max = getconservation(binalign, 'method', 'max');
if ~strcmp(cons_bkg_choice, 'rama')
    % temporarily update the background distribution
    setdbbkg('protein', []);
end

secscons = getsectors(cons(:), secopts{:}, 'refseq', alignment.refseq);
secsdiag = getsectors(diag(sca.cmat), secopts{:}, 'refseq', alignment.refseq);
disp(['Conservation calculations took ' num2str(toc(timer)) ' seconds.']);

% compare sector from SCA to 'sectors' obtained from conservation / diag(SCA)
disp('Comparison of SCA sector with conservation-based "sector":');
disp(comparesectors(secs.secs, secscons.secs));

disp('Comparison of SCA sector with diag(SCA)-based "sector":');
disp(comparesectors(secs.secs, secsdiag.secs));

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

%% Compare top eigenvector of SCA to single-site statistics

if make_figures
    figure;
    
    set(gcf, 'position', [100 200 1200 300]);
    
    % show the correlation between the top eigenvector and conservation
    subplot(131);
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
    ylabel('Top eigenvector of SCA matrix');
    beautifygraph;
    title('Conservation vs. top eigenvector of SCA matrix');
    
    % show the correlation between the top eigenvector and diag(SCA)
    subplot(132);
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
    ylabel('Top eigenvector of SCA matrix');
    beautifygraph;
    title('Diagonal vs. top eigenvector of SCA matrix');
    
    % show the correlation between diag(SCA) and conservation
    subplot(133);
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
    ylabel('Square root of diagonal of SCA matrix');
    beautifygraph;
    title('Conservation vs. diagonal of SCA matrix');
    
    set(gcf, 'name', 'Conservation, diagonal of SCA matrix, and its top eigenvector');
    preparegraph;
end

%% Check sector connectivity

[conn_secs, conninfo_secs] = checkconnectivity(secs.secs{1}, contactmap, 'pvalue', true);
[conn_cons, conninfo_cons] = checkconnectivity(secscons.secs{1}, contactmap, 'pvalue', true);
[conn_diag, conninfo_diag] = checkconnectivity(secsdiag.secs{1}, contactmap, 'pvalue', true);

disp('Connectivity of SCA sector:');
disp([num2str(conn_secs, 2) ' (p = ' num2str(conninfo_secs.p, 2) ')']);
disp('Connectivity of conserved residues:');
disp([num2str(conn_cons, 2) ' (p = ' num2str(conninfo_cons.p, 2) ')']);
disp('Connectivity of residues with highest diag(SCA):');
disp([num2str(conn_diag, 2) ' (p = ' num2str(conninfo_diag.p, 2) ')']);

%% Load the experimental data

timer = tic;
data0 = csvload('inputs/LacITable.txt', 'delimiter', ' ', 'process', false);
value_convert = containers.Map({'+', '+-', '-+', '-'}, [0 1 2 3]);
data = arrayfun(@(i) {str2double(data0.columns{1}{i}(2:end-1)), value_convert(data0.columns{2}{i}), data0.columns{1}{i}([1 end])}, ...
    1:length(data0.columns{1}), 'uniform', false);
disp(['Loading data took ' num2str(toc(timer)) ' seconds.']);

%% Collect the information for each position

timer = tic;
positions = unique(cellfun(@(c) c{1}, data));
% counting only full effects ('-', not '-+' or '+-')
effects = zeros(size(positions));
characters = blanks(length(positions));
for i = 1:length(positions)
    crtdata = data(cellfun(@(c) c{1} == positions(i), data));
    crtchar = cellfun(@(c) c{3}(1), crtdata);
    if length(unique(crtchar)) ~= 1
        warning('script:mulchars', 'Error in LacI data.');
    end
    characters(i) = crtchar(1);
    % keep only mutations to allowed residues
    mask = cellfun(@(c) ismember(c{3}(2), data_restrict), crtdata);
    effects(i) = sum(cellfun(@(c) c{2} == 3, crtdata(mask)));
end

[idxs_aln, idxs_data] = findcommonidxs(alignment.refseq.map, positions);

% test that the sequence from the data matches the one in the alignment
if ~all(alignment.data(1, idxs_aln) == characters(idxs_data))
    warning('script:badseq', 'The focal sequence from the alignment does not match that from the dataset.');
end

disp(['Preprocessing experimental data took ' num2str(toc(timer)) ' seconds.']);

%% Compare SCA and conservation to data

if make_figures
    figure;
    
    set(gcf, 'position', [100 200 800 300]);
    
    subplot(121);
    scatterfit(ev(:, 1), effects, ...
        'shapes', '.', 'sizes', 20, 'colors', 'r', 'style', {'k--', 'linewidth', 2}, ...
        'refseq', {alignment.refseq.map, positions}, 'legend', 'c');
    tmp = axis;
    axis([0 tmp(2) 0 tmp(4)]);
    xlabel('Top eigenvector of SCA matrix');
    ylabel('Functional effect');
    beautifygraph;
    title('Top eigenvector of SCA matrix vs. functional effect');
    
    subplot(122);
    scatterfit(cons, effects, ...
        'shapes', '.', 'sizes', 20, 'colors', 'r', 'style', {'k--', 'linewidth', 2}, ...
        'refseq', {alignment.refseq.map, positions}, 'legend', 'c');
    tmp = axis;
    axis([0 tmp(2) 0 tmp(4)]);
    xlabel('Conservation');
    ylabel('Functional effect');
    beautifygraph;
    title('Conservation vs. functional effect');
    
    set(gcf, 'name', 'Conservation and SCA vs. mutational effect');
    preparegraph;
end

%% Some distribution tests

disp(' ');
disp('Histogram comparisons between SCA and conservation-based sectors and data');
[idxs_aln, idxs_data] = findcommonidxs(alignment.refseq.map, positions);

hvals_all = effects(idxs_data);
hvals_sec = hvals_all(ismember(idxs_aln, secs.secs{1}));
hvals_cons = hvals_all(ismember(idxs_aln, secscons.secs{1}));

[~, ksp] = kstest2(hvals_all, hvals_sec, 0.05, 'larger');
mwwp = ranksum(hvals_all, hvals_sec, 'tail', 'left');
disp('Sector vs. all positions:');
disp(['Kolmogorov-Smirnov p = ' num2str(ksp, 2) ', Wilcoxon ranksum p = ' num2str(mwwp, 2)]);

if make_figures
    textopts = {'FontName', 'Helvetica', 'FontSize', 12};
    figure;
    subplot(2, 1, 1);
    multihist({hvals_all, hvals_sec}, linspace(0, 13, 14), ...
        'mode', 'overlap', ...
        'colors', [[0.7 0.7 0.7] ; [0 0 0]], ...
        'width', 0.9);
    tmp = axis;
    axis([-1 14 tmp(3) tmp(4)]);
    ylabel('counts');
    xlabel('mutational effect');
    title(['Sector positions (' int2str(length(hvals_sec)) ' residues)']);
    text(2.15, 125, ['Mann-Whitney p = ' num2str(mwwp, 2)], textopts{:});
    beautifygraph;
    set(gca, 'xgrid', 'on', 'ygrid', 'on', ...
        'xminorgrid', 'off', 'yminorgrid', 'off');
end

[~, ksp] = kstest2(hvals_all, hvals_cons, 0.05, 'larger');
mwwp = ranksum(hvals_all, hvals_cons, 'tail', 'left');
disp('Conservation vs. all positions:');
disp(['Kolmogorov-Smirnov p = ' num2str(ksp, 2) ', Wilcoxon ranksum p = ' num2str(mwwp, 2)]);
   
if make_figures
    subplot(2, 1, 2);
    multihist({hvals_all, hvals_cons}, linspace(0, 13, 14), ...
        'mode', 'overlap', ...
        'colors', [[0.7 0.7 0.7] ; [0 0 0]], ...
        'width', 0.9);
    tmp = axis;
    axis([-1 14 tmp(3) tmp(4)]);
    ylabel('counts');
    xlabel('mutational effect');
    title(['Top conserved positions (' int2str(length(hvals_cons)) ' residues)']);
    text(2.15, 125, ['Mann-Whitney p = ' num2str(mwwp, 2)], textopts{:});
    beautifygraph;
    set(gca, 'xgrid', 'on', 'ygrid', 'on', ...
        'xminorgrid', 'off', 'yminorgrid', 'off');
    
    set(gcf, 'units', 'normalized', 'position', [0.35 0.1 0.3 1.6/3], ...
    'name', 'Comparing sectors to conservation using histograms');
    preparegraph;
end

[~, ksp] = kstest2(hvals_sec, hvals_cons, 0.05);
mwwp = ranksum(hvals_sec, hvals_cons);
disp('Sector vs. conservation:');
disp(['Kolmogorov-Smirnov p = ' num2str(ksp, 2) ', Wilcoxon ranksum p = ' num2str(mwwp, 2)]);

%% Comparisons using Fisher exact tests

% set up some useful functions

% find set of "functionally significant" positions
fctidxs = positions(effects >= fctthreshold);
% total number of positions that appear in the experimental data
ntotpos = length(effects);

% keep only residues that also have experimental data
keepresidueswithexp = @(sec) intersect(sec, positions);
% helper function
tablefromthree = @(a, b, c, n) [a, b ; c, n - a - b - c];
generatetable0 = @(secwithexp) tablefromthree(...
    length(intersect(secwithexp, fctidxs)), length(setdiff(fctidxs, secwithexp)), ...
    length(setdiff(secwithexp, fctidxs)), ntotpos);
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
disp('Contingency tables comparing SCA or conservation-based sectors to data');

% see how well the sector does
disp('Sector:');
secftable = generatetable(secs.secs_mapped{1}{1});
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

disp(' ');
disp('Chi^2 p-value for comparing sector to conservation contingency tables:');
disp(['p = ' num2str(chi2gof2([secftable(:) consftable(:)]), 2)]);

disp(' ');
disp('Chi^2 p-value for comparing sector to diag(SCA) contingency tables:');
disp(['p = ' num2str(chi2gof2([secftable(:) diagftable(:)]), 2)]);

%% Comparisons for a range of sector sizes

[idxs_aln, idxs_data] = findcommonidxs(alignment.refseq.map, positions);

hvals_all = effects(idxs_data);

nmin = 1;
nmax = 150;

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
for n = nmin:nmax
    crtsec = getnth(1:n, getnth(2, @sort, ev(:, 1), 'descend'));
    
    crtcons = getnth(1:n, getnth(2, @sort, cons, 'descend'));
    
    crteffects_sec = hvals_all(ismember(idxs_aln, crtsec));
    crteffects_cons = hvals_all(ismember(idxs_aln, crtcons));
   
    i = n - nmin + 1;
    [~, comparison_ks(i)] = kstest2(crteffects_sec, crteffects_cons, 0.05);
    comparison_wwp(i) = ranksum(crteffects_sec, crteffects_cons);
    
    [~, sec_ks(i)] = kstest2(hvals_all, crteffects_sec, 0.05, 'larger');
    [~, cons_ks(i)] = kstest2(hvals_all, crteffects_cons, 0.05, 'larger');
    
    sec_wwp(i) = ranksum(hvals_all, crteffects_sec, 'tail', 'left');
    cons_wwp(i) = ranksum(hvals_all, crteffects_cons, 'tail', 'left');
    
    secftable = generatetable(alignment.refseq.map(crtsec));
    consftable = generatetable(alignment.refseq.map(crtcons));
    sec_pval(i) = fishertest(secftable, 'tails', 1, 'type', 'abs');
    cons_pval(i) = fishertest(consftable, 'tails', 1, 'type', 'abs');
    
    chi2_pval(i) = chi2gof2([secftable(:) consftable(:)]);
end

%% make plot for Mann Whitney p-values

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