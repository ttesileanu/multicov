% This script performs a brief analysis of SCA results on the serine
% protease alignments, as described in the supplementary information part
% of the paper. This is also the script used to generate the graphs.
%
% There are some options that can be set by defining some variables before
% running section I.0 of the script. Leave these variables undefined to get
% the default behavior.
%   'aln_choice'
%       Alignment to use. This can be
%           'pfam':     use Pfam alignment
%           'hhblits':  use HHblits alignment
%           'halabi':   use the alignment from Halabi et al.
%       (default: 'halabi')
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
%   'filtering'
%       Whether to filter the alignment. This implies first capping the
%       fraction of gaps per sequence to 20%, then eliminating sequences
%       that are more than 80% identical, and finally removing columns of
%       the alignment that contain more than 20% gaps.
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

% Tiberiu Tesileanu (2014)

%% Setup

% handle some options
if ~exist('align_choice', 'var')
    align_choice = 'halabi';
end
if ~exist('make_figures', 'var')
    make_figures = true;
end
if ~exist('filtering', 'var')
    filtering = false;
end
if ~exist('make_figures', 'var')
    make_figures = true;
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

% the positions in the experimental data are in PDB (3TGI chain E)
% coordinates, which in this case are different from Uniprot! need to do
% the translation
pdb = pdbread('pdb/3TGI.pdb');
pdb_chain = 'E';

% get the Uniprot sequence to use as a reference for numbering
[up_seq, up_info] = getdbseq({pdb, pdb_chain});

timer = tic;
        
% load the 'raw' alignment
alignment = loadfasta(['inputs/trypsin_' align_choice '.fasta'], 'protein');

% we might need a copy of the original alignment
alignment0 = alignment;

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
% update database information
alignment.refseq(1).seqdb = 'uniprot';
alignment.refseq(1).seqid = up_info.accession;
disp(['Loading alignment and truncating to PDB positions took ' num2str(toc(timer)) ' seconds.']);

if filtering
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

%% Load the experimental data

mutdata0 = csvload('inputs/halabi_mutations.csv');
mutdata = mutdata0.data;

% only interested in single mutants
mask = cellfun(@isempty, mutdata.mut2);
colnames = fieldnames(mutdata);
for i = 1:length(colnames)
    crtcol = mutdata.(colnames{i});
    mutdata.(colnames{i}) = crtcol(mask);
end

[pdbseq, pdbnames] = pdbgetseq(pdb, pdb_chain);
pdb_uniprot_map = findseqmap(pdbseq, up_seq, 'protein');
pdb_to_uniprot = containers.Map(pdbnames, pdb_uniprot_map);

% check that the wildtype amino acids mentioned in the data match the
% sequence found in the alignment
for i = 1:length(mutdata.mut)
    if strcmp(mutdata.mut{i}, 'wt')
        continue;
    end
    aa = mutdata.mut{i}(1);
    pdb_pos = mutdata.mut{i}(2:end);
    uniprot_pos = pdb_to_uniprot(pdb_pos);
    aln_pos = find(alignment.refseq.map == uniprot_pos);
    if ~isempty(aln_pos)
        if aa ~= alignment.data(1, aln_pos)
            error('script:badseq', 'Residues in experimental data do not match those in alignment sequence.');
        end
    end
end

% we prefer position indices in Uniprot instead of PDB coordinates
mutdata.pos = zeros(length(mutdata.mut), 1);

for i = 1:length(mutdata.mut)
    if ~strcmp(mutdata.mut{i}, 'wt')
        mutdata.pos(i) = pdb_to_uniprot(mutdata.mut{i}(2:end));
    end
end

%% Run SCA

% the sectors from Halabi et al., in PDB coordinates
halabisecs_pdb = ...
    {{'17', '161', '172', '176', '177', '180', '183', '187', '188', '189', '191', '192', '213', '215', '216', '220', '221', '226', '227', '228'}, ...
    {'21', '26', '46', '52', '68', '69', '71', '77', '80', '81', '104', '105', '108', '118', '123', '124', '136', '153', '157', '201', '210', '229', '237', '242', '245'}, ...
    {'19', '33', '42', '43', '55', '56', '57', '58', '102', '141', '142', '184', '194', '195', '196', '197', '198', '199', '213', '214', '216', '225'}};
% convert to Uniprot
halabisecs_up = cellfun(@(sec) cellfun(@(pos) pdb_to_uniprot(pos), sec), halabisecs_pdb, 'uniform', false);

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

%% Find the linear combinations that best approximate the Halabi sectors

% which eigenvectors to use
kmax = 1:4;

% turn sectors into masks (true = belongs to sector; false = doesn't)
halabisecs_mask = zeros(size(alignment.data, 2), length(halabisecs_up));
% Halabi et al. sectors restricted to alignment
halabisecs = cell(1, length(halabisecs_up));
for i = 1:length(halabisecs_up)
    crtsec_up = halabisecs_up{i};
    halabisecs_mask(:, i) = ismember(alignment.refseq.map, crtsec_up);
    halabisecs{i} = find(halabisecs_mask(:, i));
end

% some helper functions for displaying the plot labels below
signchars = '-+';
choosesgn = @(x) signchars(1 + ceil((sign(x) + 1)/2));
stars = {'*', ''};
starchoose = @(name) stars{isempty(name) + 1};
lincombtostr = @(fdata, names) [num2str(fdata(1), 3) '*' names{1} ' ' ...
    cell2mat(arrayfun(@(i) [choosesgn(fdata(i)) ' ' num2str(abs(fdata(i)), 3) starchoose(names{i}) names{i} ' '], 2:length(fdata), 'uniform', false))];

predictors = [ev(:, kmax) ones(size(ev, 1), 1)];
predictornames = [arrayfun(@(i) ['v' int2str(i)], kmax, 'uniform', false) {''}];
halabifits = cell(1, length(halabisecs_up));
for i = 1:length(halabifits)
    halabifits{i} = regress(halabisecs_mask(:, i), predictors);
end

%% make figures

if make_figures
    figure;
    set(gcf, 'position', [100 200 400*length(halabifits) 300]);
    
    for i = 1:length(halabifits)
        subplot(1, length(halabifits), i);
        [~, ~, tmph] = scatterfit(...
            predictors*halabifits{i}, halabisecs_mask(:, i), ...
            'shapes', '.', 'colors', 'r', ...
            'sizes', 20, 'legend', 'c', ...
            'corrtext', 'Correlation c = ', ...
            'legendloc', 'north', ...
            'showfit', true, 'showci', false, ...
            'line', [1 0], 'style', {'--k'});
        xlabel(lincombtostr(halabifits{i}, predictornames));
        ylabel(['Halabi et al. sector ' int2str(i)]);
        beautifygraph;
        set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);
        
        beautifygraph;
    end
    
    preparegraph;
end

%% get sectors from linear combinations

lincomb_secs = cell(1, length(halabifits));
for i = 1:length(lincomb_secs)
    crtlincomb = predictors*halabifits{i};
    [~, sortidxs] = sort(crtlincomb, 'descend');
    lincomb_secs{i} = sortidxs(1:length(halabisecs{i}));
end

% how similar are the sectors?
[secsim, secmap, secmat] = comparesectors(lincomb_secs, halabisecs);
disp(' ');
disp('Similarity between sectors determined from our linear combinations and those from Halabi et al.:');
for i = 1:length(lincomb_secs)
    disp([lincombtostr(halabifits{i}, predictornames) ': ' num2str(secsim(i), 2) ' similar to Halabi sector ' int2str(secmap(i))]);
end

%% Try to match experimental data with linear combinations

[idxs0_aln, idxs0_exp] = findcommonidxs(alignment.refseq.map, mutdata.pos);

% get those portions of the eigenvectors that match the available
% mutagenesis data
ev_masked = ev(idxs0_aln, :);
% get those portions of the data that match the alignment
mutdata_masked = mutdata;
names = fieldnames(mutdata);
for i = 1:length(names)
    crtcol = mutdata.(names{i});
    mutdata_masked.(names{i}) = crtcol(idxs0_exp);
end

% fit!
predictors_masked = [ev_masked(:, kmax) ones(size(ev_masked, 1), 1)];
fit_tm = regress(mutdata_masked.tm, predictors_masked);
fit_kcatkm = regress(log10(mutdata_masked.kcatkm), predictors_masked);

%% make figures

if make_figures
    figure;
    [~, ~, tmph] = scatterfit(...
        predictors_masked*fit_tm, mutdata_masked.tm, ...
        'shapes', '.', 'colors', 'r', ...
        'sizes', 20, 'legend', 'c', ...
        'corrtext', 'Correlation c = ', ...
        'legendloc', 'north', ...
        'showfit', true, 'showci', false, ...
        'line', [1 0], 'style', {'--k'});
    xlabel(lincombtostr(fit_tm, predictornames));
    ylabel('T_m');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);
    tmp = get(gcf, 'position');
    set(gcf, 'position', [tmp(1:2) tmp(3:4)*0.75]);
    preparegraph;
    
    figure;
    [~, ~, tmph] = scatterfit(...
        predictors_masked*fit_kcatkm, log10(mutdata_masked.kcatkm), ...
        'shapes', '.', 'colors', 'r', ...
        'sizes', 20, 'legend', 'c', ...
        'corrtext', 'Correlation c = ', ...
        'legendloc', 'north', ...
        'showfit', true, 'showci', false, ...
        'line', [1 0], 'style', {'--k'});
    xlabel(lincombtostr(fit_kcatkm, predictornames));
    ylabel('log_{10}(k_{cat}/k_m)');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);
    tmp = get(gcf, 'position');
    set(gcf, 'position', [tmp(1:2) tmp(3:4)*0.75]);
    preparegraph;
end

%% How does conservation fare with the experimental data?

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

cons_masked = cons(idxs0_aln);

%% make figures

if make_figures
    figure;
    set(gcf, 'position', [100 200 800 300]);
    
    subplot(1, 2, 1);
    [~, ~, tmph] = scatterfit(...
        cons_masked, mutdata_masked.tm, ...
        'shapes', '.', 'colors', 'r', ...
        'sizes', 20, 'legend', 'c', ...
        'corrtext', 'Correlation c = ', ...
        'legendloc', 'north', ...
        'showfit', true, 'showci', false, ...
        'style', {'-k'});
    xlabel('Conservation');
    ylabel('T_m');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);

    subplot(1, 2, 2);
    [~, ~, tmph] = scatterfit(...
        cons_masked, log10(mutdata_masked.kcatkm), ...
        'shapes', '.', 'colors', 'r', ...
        'sizes', 20, 'legend', 'c', ...
        'corrtext', 'Correlation c = ', ...
        'legendloc', 'north', ...
        'showfit', true, 'showci', false, ...
        'style', {'-k'});
    xlabel('Conservation');
    ylabel('log_{10}(k_{cat}/k_m)');
    beautifygraph;
    set(tmph.legend, 'fontname', 'avantgarde', 'fontsize', 12);

    preparegraph;
end