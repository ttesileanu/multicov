% This script performs the analysis of the SCA method on an alignment of
% voltage-sensitive domains of potassium channels, as described in the
% paper. This is also the script used to generate the graphs.
%
% The basic purpose of this script is to show how the sector positions
% and the top conserved positions compare to the mutagenesis data from
% Li-Smerin et al. 2000 and Lee et al. 2009.
%
% There are some options that can be set by defining some variables before
% running the first section of the script. Leave these variables undefined
% to get the default behavior.
%   'aln_choice'
%       Which alignment to use.
%         'lee2009': alignment from Lee et al. 2009
%         'hhblits_e10': HHblits alignment, e-value: e10
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
%   'dataset'
%       Which dataset to use. Options: 'lismerin', 'lee2009'.
%       (default: 'lismerin')
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

if ~exist('dataset', 'var')
    dataset = 'lismerin';
end
if ~exist('sca_type', 'var')
    sca_type = 'projection';
end
if ~exist('sec_fraction', 'var')
    sec_fraction = 0.25;
end
if ~exist('filtering', 'var')
    filtering = false;
end
if ~exist('make_figures', 'var')
    make_figures = true;
end
if ~exist('aln_choice', 'var')
    aln_choice = 'hhblits_e10';
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
% Uniprot IDs for the two proteins used in Li-Smerin et al. and in Lee et al.
up_id_kv12 = 'KCNA2_RAT';
up_id_kv21 = 'KCNB1_RAT';

[kv12_seq, kv12_info] = getdbseq('uniprot', up_id_kv12);
[kv21_seq, kv21_info] = getdbseq('uniprot', up_id_kv21);
% for some reason, need a shift to match with Li-Smerin et al.'s coordinates
kv21_seq = kv21_seq(5:end);

switch aln_choice
    case 'lee2009'
        % don't ignore "insert" states -- if we want to filter out gappy columns,
        % we can do it later
        alignment = loadfasta('inputs/kv12alignment_lee2009.fasta', 'protein', ...
            'keepmask', @(v) true(size(v)), 'postprocess', @(v) upper(v));
        
        % keep only the range of residues that we have data for
        datarange = 162:314;
        
        % find the focal sequence and put it in first place
        alignment0 = alignment;
        alignment = seqsearch(alignment, 'seq', kv12_seq(datarange));
        alignment = seqtruncate(alignment, 'seq', kv12_seq(datarange), 'posnames', datarange);
        alignment.refseq(1).seqdb = 'uniprot';
        alignment.refseq(1).seqid = kv12_info.accession;
    case 'hhblits_e10'
        alignment = loadfasta('inputs/KCNB1_RAT_e10_n2_m40.a3m', 'protein');
        
        % keep only the range of residues that we have data for
        datarange = 185:311;
        
        % find the focal sequence and put it in first place
        alignment0 = alignment;
        alignment = seqsearch(alignment, 'annot', kv21_info.id, 'seq', kv21_seq(datarange));
        alignment = seqtruncate(alignment, 'seq', kv21_seq(datarange), 'posnames', datarange);
        alignment.refseq(1).seqdb = 'uniprot';
        alignment.refseq(1).seqid = kv21_info.accession;
end
disp(['Loading alignment and truncating to relevant range took ' num2str(toc(timer)) ' seconds.']);

% filter the alignment, if required
if filtering
    timer = tic; %#ok<*UNRCH>
    % eliminate fragments
    alignment = filterrows(alignment, 0.2);
    % eliminate near-duplicates
    alignment = eliminatesimilarseqs(alignment, 0.8);
    % eliminate positions that most sequences don't have
    alignment = filtercolumns(alignment, 0.2);
    
%    save(['generated/kv_' aln_choice '_filtered.mat'], 'alignment');
    disp(['Filtering the alignment took ' num2str(toc(timer)) ' seconds.']);
end

% decide on a sector size as a number of residues
if sec_fraction < 1
    sec_count = round(size(alignment.data, 2)*sec_fraction);
else
    sec_count = sec_fraction;
end

%% What does SCA make of this?

if ~strcmp(sca_bkg_choice, 'rama')
    % temporarily update the background distribution
    setdbbkg('protein', sca_bkg_choice);
end
sca = getsca(alignment, 'method', sca_type);
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

%% Get conservation

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

secscons = getsectors(cons(:), secopts{:}, 'refseq', alignment.refseq);
secsdiag = getsectors(diag(sca.cmat), secopts{:}, 'refseq', alignment.refseq);
disp(['Conservation calculations took ' num2str(toc(timer)) ' seconds.']);

% compare SCA sector to sectors obtained from conservation / diag(SCA)
disp('Comparison of SCA sector with conservation-based "sector":');
disp(comparesectors(secs.secs, secscons.secs));

disp('Comparison of SCA sector with diag(SCA)-based "sector":');
disp(comparesectors(secs.secs, secsdiag.secs));

%% Compare different definitions for conservation

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

%% Load the experimental data

timer = tic;
data0_lismerin = csvload('inputs/smerinetal_alaninescan.csv');
data0_lee2009 = csvload('inputs/lee2009_data.csv');
disp(['Loading data took ' num2str(toc(timer)) ' seconds.']);

% NOTE:
% 1) compared to the table in Li-Smerin et al., there is one important
% difference: the \Delta\Delta G for P248 seems to have been calculated or
% transcribed wrongly in the paper (the \Delta G value is correctly equal
% to the 0.2389*z*F*V_{50}, as described, but \Delta\Delta G does not
% follow from that); here we calculate \Delta\Delta G in the Matlab script,
% so this issue does not appear
% 2) because we recalculate \Delta\Delta G using the truncated numbers
% found in the Li-Smerin et al. paper, it is not surprising if some
% uncertainties in the calculations lead to different values of
% \Delta\Delta G; the one noted above is too large to be accounted for by
% this, but one case in which uncertainties in calculations seem to be
% important is V292, where Li-Smerin gives \Delta\Delta G = 1.5 kcal/mol,
% while our calculation yields 1.6. The values are, however, within one
% standard deviation (equal to 0.1 in this case).

%% Map Li-Smerin et al. and Lee et al. data to alignment numbering

if strcmp(aln_choice, 'lee2009')
    % find mapping between Kv2.1 numbering and Kv1.2 numbering
    [~, kv21_idx] = seqsearch(alignment, 'seq', kv21_seq, 'swaptofirst', false, 'verbose', false);
    kv21_map = zeros(1, size(alignment.data, 2));
    % each element in this vector maps an alignment position to a Kv2.1 index
    kv21_map(alignment.data(kv21_idx, :) ~= '-') = findseqmap(alignment.data(kv21_idx, :), kv21_seq, 'protein');
else
    % we still need kv21_map for what is happenning below, but since the
    % hhblits alignment is focused on Kv2.1, this is trivial
    kv21_map = alignment.refseq.map;
    kv21_idx = 1;
end
    
lismerin_effects = nan(size(kv21_map));
v50_effects = nan(size(kv21_map));
z_effects = nan(size(kv21_map));
lismerin_std = nan(size(kv21_map));
lismerin_idxs = cellfun(@(s) str2double(s(2:end)), data0_lismerin.data.mut(2:end));
%lismerin_mask = (abs(data0_lismerin.data.dg(2:end)) >= 1);

if strcmp(aln_choice, 'lee2009')
    lee2009_effects = nan(size(kv21_map));
    lee2009_effects(data0_lee2009.data.pos_kv12 - datarange(1) + 1) = data0_lee2009.data.effect_ala;
    datamask_lee2009 = ~isnan(lee2009_effects);
    lee2009_effects_nonan = lee2009_effects;
    lee2009_effects_nonan(~datamask_lee2009) = -1;
end

for i = 1:length(lismerin_idxs)
    idx = lismerin_idxs(i);
    revidx = find(kv21_map == idx);
    if ~isempty(revidx)
        chr_aln = alignment.data(kv21_idx, revidx);
        mut_aln = data0_lismerin.data.mut{1 + i}(1);
        if chr_aln ~= mut_aln
            warning(['Expected ' mut_aln int2str(idx) ', found ' chr_aln ' at alignment position ' int2str(revidx) '.']);
        end
        lismerin_effects(revidx) = data0_lismerin.data.dg(1 + i) - data0_lismerin.data.dg(1);
        v50_effects(revidx) = data0_lismerin.data.v50(1 + i);
        z_effects(revidx) = data0_lismerin.data.z(1 + i);
        lismerin_std(revidx) = data0_lismerin.data.std_ddg(1 + i);
    end
end

datamask_lismerin = ~isnan(lismerin_effects);

lismerin_effects_bin = nan(size(lismerin_effects));
lismerin_effects_bin(datamask_lismerin) = (abs(lismerin_effects(datamask_lismerin)) >= 1 - 1e-6);

lismerin_effects_nonan = lismerin_effects_bin;
lismerin_effects_nonan(~datamask_lismerin) = -1;


% Does the data from the two papers match? (They should; Lee et al. use
% Li-Smerin et al. as a souce of experimental data.)

% Kind of. Exceptions are the following:
% 1) There are many places where Li-Smerin et al. show data, but this is
% not shown in the Lee paper.
% 2) There are also two positions where the Lee paper shows Ala data, but
% it's not present in the Li-Smerin table.
% 3) Lee 2009 seems to have used the \Delta\Delta G from the Li-Smerin et
% al. table, and this is wrong for P248.
% 4) For I294, \Delta\Delta G is -0.9, but it was counted as significant by
% Lee 2009. It is not clear why.

%% choose a dataset

switch dataset
    case 'lismerin'
        effects = abs(lismerin_effects);
        effects_bin = lismerin_effects_bin;
        datamask = datamask_lismerin;
    case 'lee2009'
        if strcmp(aln_choice, 'lee2009')
            effects = lee2009_effects;
            effects_bin = effects;
            datamask = datamask_lee2009;
        else
            error('Cannot use Lee 2009 dataset with hhblits alignment.');
        end
    otherwise
        error(['Unrecognized dataset option ' dataset '.']);
end

%% Compare SCA and conservation to data

if make_figures
    figure;
    
    set(gcf, 'position', [100 200 800 300]);
    
    subplot(121);
    scatterfit(ev(datamask, 1), effects(datamask), ...
        'shapes', '.', 'sizes', 20, 'colors', 'r', 'style', {'k--', 'linewidth', 2}, ...
        'legend', 'c');
    tmp = axis;
    axis([0 tmp(2) 0 tmp(4)]);
    xlabel('Top eigenvector of SCA matrix');
    ylabel('Functional effect');
    beautifygraph;
    title('Top eigenvector of SCA matrix vs. functional effect');
    
    subplot(122);
    scatterfit(cons(datamask), effects(datamask), ...
        'shapes', '.', 'sizes', 20, 'colors', 'r', 'style', {'k--', 'linewidth', 2}, ...
        'legend', 'c');
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

makebinvec = @(pos, n) full(sparse(pos, ones(1, length(pos)), ones(1, length(pos)), n, 1));
secs1_mask = makebinvec(secs.secs{1}, size(alignment.data, 2));
secscons1_mask = makebinvec(secscons.secs{1}, size(alignment.data, 2));

effects_sec = effects(datamask(:) & secs1_mask);
effects_cons = effects(datamask(:) & secscons1_mask);

[~, ksp] = kstest2(effects(datamask), effects_sec, 0.05, 'larger');
mwwp = ranksum(effects(datamask), effects_sec, 'tail', 'left');
disp('Sector vs. all positions:');
disp(['Kolmogorov-Smirnov p = ' num2str(ksp, 2) ', Wilcoxon ranksum p = ' num2str(mwwp, 2)]);

if make_figures
    textopts = {'FontName', 'Helvetica', 'FontSize', 12};
    figure;
    subplot(2, 1, 1);
    multihist({effects(datamask), effects_sec}, linspace(0, 13, 14), ...
        'mode', 'overlap', ...
        'colors', [[0.7 0.7 0.7] ; [0 0 0]], ...
        'width', 0.9);
    tmp = axis;
    axis([-1 14 tmp(3) tmp(4)]);
    ylabel('counts');
    xlabel('mutational effect');
    title(['Sector positions (' int2str(length(effects_sec)) ' residues)']);
    text(2.15, 0.25*tmp(3) + 0.75*tmp(4), ['Mann-Whitney p = ' num2str(mwwp, 2)], ...
        textopts{:});
    beautifygraph;
    set(gca, 'xgrid', 'on', 'ygrid', 'on', ...
        'xminorgrid', 'off', 'yminorgrid', 'off');
end

[~, ksp] = kstest2(effects(datamask), effects_cons, 0.05, 'larger');
mwwp = ranksum(effects(datamask), effects_cons, 'tail', 'left');
disp('Conservation vs. all positions:');
disp(['Kolmogorov-Smirnov p = ' num2str(ksp, 2) ', Wilcoxon ranksum p = ' num2str(mwwp, 2)]);

if make_figures
    subplot(2, 1, 2);
    multihist({effects(datamask), effects_cons}, linspace(0, 13, 14), ...
        'mode', 'overlap', ...
        'colors', [[0.7 0.7 0.7] ; [0 0 0]], ...
        'width', 0.9);
    tmp = axis;
    axis([-1 14 tmp(3) tmp(4)]);
    ylabel('counts');
    xlabel('mutational effect');
    title(['Top conserved positions (' int2str(length(effects_sec)) ' residues)']);
    text(2.15, 0.25*tmp(3) + 0.75*tmp(4), ['Mann-Whitney p = ' num2str(mwwp, 2)], ...
        textopts{:});
    beautifygraph;
    set(gca, 'xgrid', 'on', 'ygrid', 'on', ...
        'xminorgrid', 'off', 'yminorgrid', 'off');
    
    set(gcf, 'units', 'normalized', 'position', [0.35 0.1 0.3 1.6/3], ...
        'name', 'Comparing sectors to conservation using histograms');
    preparegraph;
end

[~, ksp] = kstest2(effects_sec, effects_cons, 0.05);
mwwp = ranksum(effects_sec, effects_cons);
disp('Sector vs. conservation:');
disp(['Kolmogorov-Smirnov p = ' num2str(ksp, 2) ', Wilcoxon ranksum p = ' num2str(mwwp, 2)]);

%% Comparison using Fisher's exact tests

% set up some useful functions

% find set of "functionally significant" positions
% ugh, nan counts as different from zero...
fctidxs = find(effects_bin ~= 0 & datamask);
% total number of positions that appear in the experimental data
ntotpos = sum(datamask);

% keep only residues that also have experimental data
keepresidueswithexp = @(sec) sec(datamask(sec));
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
secftable = generatetable(secs.secs{1});
disp(generateformattedtab(secftable));
disp(generatepvaltext(secftable));
disp(' ');

% make a Fisher's exact table/test for conservation
disp('Conservation:');
consftable = generatetable(secscons.secs{1});
disp(generateformattedtab(consftable));
disp(generatepvaltext(consftable));
disp(' ');
        
% make a Fisher's exact table/test for diag(SCA)
disp('Diagonal of SCA matrix:');
diagftable = generatetable(secsdiag.secs{1});
disp(generateformattedtab(diagftable));
disp(generatepvaltext(diagftable));

disp(' ');
disp('Chi^2 p-value for comparing sector to conservation contingency tables:');
disp(['p = ' num2str(chi2gof2([secftable(:) consftable(:)]), 2)]);

disp(' ');
disp('Chi^2 p-value for comparing sector to diag(SCA) contingency tables:');
disp(['p = ' num2str(chi2gof2([secftable(:) diagftable(:)]), 2)]);

%% Comparison for a range of sector sizes

nmin = 1;
nmax = 75;

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
    crtsec_mask = makebinvec(crtsec, size(alignment.data, 2));
    
    crtcons = getnth(1:n, getnth(2, @sort, cons, 'descend'));
    crtcons_mask = makebinvec(crtcons, size(alignment.data, 2));
    
    crteffects_sec = effects(datamask(:) & crtsec_mask);
    crteffects_cons = effects(datamask(:) & crtcons_mask);
   
    i = n - nmin + 1;
    if ~isempty(crteffects_sec)
        if ~isempty(crteffects_cons)
            [~, comparison_ks(i)] = kstest2(crteffects_sec, crteffects_cons, 0.05);
            comparison_wwp(i) = ranksum(crteffects_sec, crteffects_cons);
        else
            comparison_ks(i) = nan;
            comparison_wwp(i) = nan;
        end
        [~, sec_ks(i)] = kstest2(effects(datamask), crteffects_sec, 0.05, 'larger');
        sec_wwp(i) = ranksum(effects(datamask), crteffects_sec, 'tail', 'left');    
    else
        sec_ks(i) = nan;
        sec_wwp(i) = nan;
    end
    
    if ~isempty(crteffects_cons)
        [~, cons_ks(i)] = kstest2(effects(datamask), crteffects_cons, 0.05, 'larger');
        cons_wwp(i) = ranksum(effects(datamask), crteffects_cons, 'tail', 'left');
    else
        cons_ks(i) = nan;
        cons_wwp(i) = nan;
    end
    
    
    secftable = generatetable(crtsec);
    consftable = generatetable(crtcons);
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

%% Looking for linear combinations of eigenvectors that match the data

% get those portions of the eigenvectors that match the available
% mutagenesis data
ev_masked = ev(datamask_lismerin, :);
v50_masked = v50_effects(datamask_lismerin);
z_masked = z_effects(datamask_lismerin);

% predictors are the top 3 eigenvectors of the SCA matrix
predictors = [ev_masked(:, 1:3) ones(size(ev_masked, 1), 1)];
[fit_v50, int_v50, ~, ~, stats_v50] = regress(v50_masked(:), predictors);
[fit_z, int_z, ~, ~, stats_z] = regress(z_masked(:), predictors);

% some functions useful for displaying the labels
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
        predictors*fit_v50, v50_masked, ...
        'shapes', '.', 'colors', 'r', ...
        'sizes', 20, 'legend', 'c', ...
        'corrtext', 'Correlation c = ', ...
        'legendloc', 'north', ...
        'showfit', true, 'showci', false, ...
        'line', [1 0], 'style', {'--k'});
    xlabel(lincombtostr(fit_v50, {'v1', 'v2', 'v3', ''}));
    ylabel('V_{50}');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);
    tmp = get(gcf, 'position');
    set(gcf, 'position', [tmp(1:2) tmp(3:4)*0.75]);
    preparegraph;
    
    figure;
    [~, ~, tmph] = scatterfit(...
        predictors*fit_z, z_masked, ...
        'shapes', '.', 'colors', 'r', ...
        'sizes', 20, 'legend', 'c', ...
        'corrtext', 'Correlation c = ', ...
        'legendloc', 'north', ...
        'showfit', true, 'showci', false, ...
        'line', [1 0], 'style', {'--k'});
    xlabel(lincombtostr(fit_z, {'v1', 'v2', 'v3', ''}));
    ylabel('z');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);
    tmp = get(gcf, 'position');
    set(gcf, 'position', [tmp(1:2) tmp(3:4)*0.75]);
    preparegraph;
end