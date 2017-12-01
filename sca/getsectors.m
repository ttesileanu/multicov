function [secs, extra] = getsectors(comps, varargin)
% GETSECTORS Identify statistical sectors based on spectral analysis.
%   secdata = GETSECTORS(comps) identifies sectors given a set of vectors
%   (e.g., eigenvectors of the SCA matrix). The number of sectors is given
%   by the number of columns in comps. The return value 'secdata' is a
%   structure with fields
%     'type'  -- sectors
%     'secs'  -- the sector definition, as a cell array of vectors of
%                indices;
%     'comps' -- the component matrix used to get the sectors.
%
%   When reference sequence mapping is enabled (see 'refseq' option),
%   secdata also contains a field 'secs_mapped'. This is a cell array with
%   as many entries as there are fields in the reference sequence, and each
%   of these entries is a cell array of vectors, one vector of mapped
%   residues for each sector.
%
%   GETSECTORS({stats, n}) uses the top n eigenvectors of the 'cmat' field
%   of the statistics structure 'stats'. If n is a vector instead of an
%   integer, it is used as a set of indices identifying which eigenvectors
%   of the 'cmat' matrix are to be used.
%
%   There are a number of options allowed for this function. As a general
%   rule, most of the options can either contain one choice that applies to
%   all the sectors, or an array of choices, one for each sector.
%
%   Options:
%     'cutoff' <x>/<v>
%       Cutoff to use for finding sector positions. The meaning of the
%       cutoff depends on the '-cutoff-type', and on the number of '-tails'.
%       See those options for details.
%       (default: depends on 'cutofftype':
%           0.05 for 'magnitude',
%           2 for 'range',
%           0.9 for 'cdf',
%           0.05 for 'pdf',
%           20 for 'count')
%     'distrib' <s>/<c>
%       Probability distribution to fit to histogram of component values.
%       Use 'none' for no fit -- use sample statistics or histograms.
%       (default: 'none')
%     'figures' <b>
%       Whether to make histograms and plots showing the way in which the
%       sectors were selected. A scatter plot is generated between the
%       components of vectors 1 and 2, one between vectors 3 and 4, etc. If
%       the number of vectors is odd, then there will be a scatter plot
%       between the first and last vectors. Note that in the scatter plots,
%       the sectors are drawn in order, so residues that belong to e.g. the
%       last sector will be colored in the last sector's color even if they
%       also belong to other sectors.
%       (default: false)
%     'fixmu' <x>/<v>
%       The thresholds when 'method' is 'range' are imposed with respect to
%       the mean of the distribution of components. Using the 'fixmu' option,
%       one can instead force a different reference point to be used
%       (effectively setting the value of the mean). Note that if 'distrib'
%       is 'none', the standard deviations will be calculated around this
%       forced value of the mean, but this currently does not work for
%       other values of 'distrib'. Empty array means infer from data.
%       (default: [])
%     'fixstd' <x>/<v>
%       Same as above, but fix the standard deviation value used when
%       'method' is 'range'. Again, an empty array means infer from data.
%       (default: [])
%     'method' <s>/<c>
%       What kind of threshold to put on sector components. This can be
%        'magnitude'
%           Put the threshold directly on the magnitude of the components,
%           either in absolute value (if 'tails' is 2), or on the signed
%           value (if 'tails' is 1).
%        'range'
%           Components at least 'cutoff' standard deviations away from the
%           mean are considered part of the sector. If this is one-tailed,
%           only deviations in the direction of the tail are considered
%           significant. Otherwise, the selection is done symmetrically
%           around the mean.
%        'cdf'
%           The threshold is chosen so that a fraction equal to 'cutoff' of
%           the distribution of component values is non-significant, i.e.,
%           does not belong to the sector. If this is one-tailed, the CDF
%           cut is placed in the direction of the tail. Otherwise, the cut
%           is placed symmetrically on the two sides.
%        'pdf'
%           The components with a probability smaller than 'cutoff' are
%           considered sector residues. This only works if 'distrib' is not
%           'none'. If 'tails' is 1, only residues in the tail are
%           considered significant. Otherwise, the same probability
%           threshold is applied symmetrically to the two tails.
%        'count'
%           Make sectors of exactly the size given by 'cutoff'.
%       Note that CDFs, means, and standard deviations are either calculated
%       exactly from the fitted distributions (if 'distrib' is differet from
%       'none'), or otherwise are evaluated directly from the components.
%       (default: 'range')
%     'overlaps' <b>
%       If this is true, overlaps are allowed between the sectors. If it is
%       false, a residue is considered part of a sector only if it does not
%       fulfill the thresholds for any other component.
%       (default: true)
%     'refseq' <struct>
%       If the function is used with the 'comps' form (i.e., missing refseq
%       information), then this option can be used to provide reference
%       sequence information, in the same format as is used for alignments.
%       This can also be used to override the refseq information from the
%       SCA structure. Use an empty array to disable mapping to reference
%       sequence.
%       (default: [])
%     'shade' <b>
%       If this is true, each residue's saturation in the plots made when
%       'figures' is true depends on how large the corresponding component
%       is, i.e., how 'much' the residue belongs to the sector.
%       (default: false)
%     'tails' <n>/<v>
%       Number of tails to use for thresholding -- should be 1 or 2. If 1,
%       the threshold is placed only on the positive side of comps. If the
%       form of the function is used that calculates the components from the
%       'cmat' matrix, the sign of the components is chosen using the
%       standard method in eigsorted. See details for the types of cutoffs
%       that are available under 'cutofftype'.
%       (default: 2)
%     'title' <s>
%       A string to be included in the titles of all figures created if
%       'figures' is true.
%       (default: '')
%
%   See also: GETSCA, EIGSORTED.

% Tiberiu Tesileanu (2012-2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('cutoff', [], @(v) isvector(v) && isnumeric(v));
parser.addParamValue('distrib', 'none', @(v) (ischar(v) && isvector(v)) || iscell(v));
parser.addParamValue('figures', false, @(b) isscalar(b) && islogical(b));
parser.addParamValue('fixmu', [], @(v) (isvector(v) && isnumeric(v)) || iscell(v));
parser.addParamValue('fixstd', [], @(v) (isvector(v) && isnumeric(v)) || iscell(v));
parser.addParamValue('method', 'range', @(v) (ischar(v) && isvector(v)) || iscell(v));
parser.addParamValue('overlaps', true, @(b) isscalar(b) && islogical(b));
parser.addParamValue('refseq', [], @(s) isstruct(s));
parser.addParamValue('shade', false, @(b) isscalar(b) && islogical(b));
parser.addParamValue('tails', 2, @(v) isvector(v) && isnumeric(v));
parser.addParamValue('title', '', @(s) ischar(s) && isvector(s));

% parse
parser.parse(varargin{:});
params = parser.Results;

if iscell(comps)
    % we are asked to calculate the components on our own
    sca = comps{1};
    n = comps{2};
    
    if isempty(params.refseq)
        params.refseq = sca.refseq;
    end
    
    if isscalar(n)
        myrange = 1:n;
    else
        myrange = n;
    end
    
    [ev_sca, ~] = eigsorted(sca.cmat);
    
    comps = ev_sca(:, myrange);
end

[Npos, nfit] = size(comps);

% extend distrib to one for each sector
if ~iscell(params.distrib)
    params.distrib = repmat({params.distrib}, nfit, 1);
elseif length(params.distrib) ~= nfit
    error([mfilename ':baddistvec'], 'Mismatch between number of sectors and length of ''distrib'' array.');
end
params.distrib = params.distrib(:);
% extend method to one for each sector
if ~iscell(params.method)
    params.method = repmat({params.method}, nfit, 1);
elseif length(params.method) ~= nfit
    error([mfilename ':badctvec'], 'Mismatch between number of sectors and length of ''method'' array.');
end
if ~all(ismember(params.method, {'range', 'pdf', 'cdf', 'magnitude', 'count'}))
    error([mfilename ':badct'], 'Unrecognized method.');
end
params.method = params.method(:);
% extend tails to one for each sector
if isscalar(params.tails)
    params.tails = repmat(params.tails, nfit, 1);
elseif length(params.tails) ~= nfit
    error([mfilename ':badctvec'], 'Mismatch between number of sectors and length of ''tails'' array.');
end
if ~all(params.tails == 1 | params.tails == 2)
    error([mfilename ':badtails'], 'Tails should be 1 or 2.');
end
params.tails = params.tails(:);
% check whether we've been given explicit cutoffs
if isempty(params.cutoff)
    % no, use the defaults
    params.cutoff = zeros(nfit, 1);
    for i = 1:length(params.method)
        switch params.method{i}
            case 'magnitude'
                params.cutoff(i) = 0.05;
            case 'range'
                params.cutoff(i) = 2;
            case 'cdf'
                params.cutoff(i) = 0.9;
            case 'pdf'
                params.cutoff(i) = 0.05;
            case 'count'
                params.cutoff(i) = 20;
        end
    end
else
    % extend cutoffs to one for each sector
    if isscalar(params.cutoff)
        params.cutoff = repmat(params.cutoff, nfit, 1);
    elseif length(params.cutoff) ~= nfit
        error([mfilename ':badctvec'], 'Mismatch between number of sectors and length of ''cutoff'' array.');
    end
end
params.cutoff = params.cutoff(:);
% extend fixmu for each sector
if ~iscell(params.fixmu)
    if isempty(params.fixmu) || isscalar(params.fixmu)
        params.fixmu = repmat({params.fixmu}, nfit, 1);
    else
        params.fixmu = num2cell(params.fixmu);
    end
end
if length(params.fixmu) ~= nfit
    error([mfilename ':badfixmu'], 'Mismatch between number of sectors and length of ''fixmu'' array.');
end
params.fixmu = params.fixmu(:);
% extend fixstd for each sector
if ~iscell(params.fixstd)
    if isempty(params.fixstd) || isscalar(params.fixstd)
        params.fixstd = repmat({params.fixstd}, nfit, 1);
    else
        params.fixstd = num2cell(params.fixstd);
    end
end
if length(params.fixstd) ~= nfit
    error([mfilename ':badfixstd'], 'Mismatch between number of sectors and length of ''fixstd'' array.');
end
params.fixstd = params.fixstd(:);

% prepare the figures, if required
if params.figures
    figure;
    set(gcf, 'Units', 'normalized', 'Position', [0 0.6 min(1, 0.2*nfit) 0.4], 'Name', 'Sector Definition');
end

% start searching for sectors!
secs0 = cell(nfit, 1);
cutoffs = repmat([-inf ; inf], 1, nfit);
for i = 1:nfit
    vec = comps(:, i);
    
    % draw the histogram
    if params.figures
        % Freedman-Diaconis rule
        binwidth = 2*iqr(vec)*Npos^(-1/3);
        nbins = ceil(range(vec)/binwidth);
        if nbins < 8 % don't make the bins too large
            nbins = 8;
        end

        [yhist, xhist] = hist(comps(:, i), nbins);
        subplot(1, nfit, i);
        bar(xhist, yhist, 'k');
        crttitle = ['sector ' int2str(i)];
        if ~isempty(params.title)
            crttitle = [params.title ', ' crttitle]; %#ok<AGROW>
        end
        title(crttitle, 'FontSize', 12);
        hold on;
        grid on;
    end
    
    % fit a distribution, if asked for that
    if ~strcmp(params.distrib{i}, 'none')
        fit = fitdist(vec, params.distrib{i});
        % get the PDF of the distribution
        % XXX make the number of points configurable? larger/smaller?
        xpts = linspace(min(vec), max(vec), 4096);
        mypdf = pdf(fit, xpts);
        % plot it, if required
        if params.figures
            % need some scaling to match to histogram
            scaled_pdf = (xhist(2) - xhist(1))*Npos*mypdf;
            plot(xpts, scaled_pdf, 'r-', 'LineWidth', 2);
        end
        
        % for possible later usage, get the mean and standard deviation
        % but only if they're not fixed by the user
        if isempty(params.fixmu{i})
            mu = mean(fit);
        else
            mu = params.fixmu{i};
        end
        % XXX when we fix mu, it would make sense to calculate sigma around
        % the forced mean
        if isempty(params.fixstd{i})
            sigma = std(fit);
        else
            sigma = params.fixstd{i};
        end
    else
        % for later usage, get the mean and standard deviation
        if isempty(params.fixmu{i})
            mu = mean(vec);
        else
            mu = params.fixmu{i};
        end
        if isempty(params.fixstd{i})
            % if we changed the mean, make sure to calculate the standard
            % deviation around the mean
            if isempty(params.fixmu{i})
                sigma = std(vec);
            else
                sigma = sqrt(sum((vec - mu).^2) / (length(vec) - 1));
            end
        else
            sigma = params.fixstd{i};
        end
    end
    
    % find the cutoffs
    switch params.method{i}
        case 'magnitude'
            cutoffs(2, i) = params.cutoff(i);
            if params.tails(i) > 1
                cutoffs(1, i) = -params.cutoff(i);
            end
        case 'range'
            cutoffs(2, i) = mu + params.cutoff(i)*sigma;
            if params.tails(i) > 1
                cutoffs(1, i) = mu - params.cutoff(i)*sigma;
            end
        case 'cdf'
            if ~strcmp(params.distrib{i}, 'none')
                % get the CDF
                mycdf = cdf(fit, xpts);
                if params.tails(i) == 1
                    [~, x1] = find(mycdf >= params.cutoff(i), 1);
                    if ~isempty(x1)
                        cutoffs(2, i) = xpts(x1);
                    end
                else
                    % need half the cutoff for each tail
                    cutoff12 = 1 - (1 - params.cutoff(i))/2;
                    [~, x1] = find(mycdf >= 1 - cutoff12, 1);
                    [~, x2] = find(mycdf >= cutoff12, 1);
                    if ~isempty(x1)
                        cutoffs(1, i) = xpts(x1);
                    end
                    if ~isempty(x2)
                        cutoffs(2, i) = xpts(x2);
                    end
                end
            else
                if params.tails(i) == 1
                    sorted = sort(vec);
                    cutoffs(2, i) = sorted(ceil(params.cutoff(i)*length(vec)));
                else
                    sorted = sort(vec);
                    cutoff12 = 1 - (1 - params.cutoff(i))/2;
                    cutoffs(1, i) = sorted(ceil((1 - cutoff12)*length(vec)));
                    cutoffs(2, i) = sorted(ceil(cutoff12*length(vec)));
                end
            end
        case 'pdf'
            if strcmp(params.distrib{i}, 'none')
                error([mfilename ':needdist'], 'Need a non-trivial distribution when ''method'' is ''pdf''.');
            end
            cutoffs(2, i) = xpts(find(mypdf > params.cutoff(i), 1, 'last'));
            if params.tails(i) > 1
                cutoffs(1, i) = xpts(find(mypdf > params.cutoff(i), 1, 'first'));
            end
        case 'count'
        otherwise
            error([mfilename ':internal'], 'This should never happen: bad method.');
    end
    
    if ~strcmp(params.method{i}, 'count')
        % we obtain the indices of sector positions
        crtsec = find(vec <= cutoffs(1, i) | vec >= cutoffs(2, i));
        % let's sort them in a more informative order -- based on distance from
        % cutoff for single-sided tests, or distance from mean of cutoffs for
        % two-sided
        if isinf(cutoffs(1, i))
            secs0{i} = crtsec(getnth(2, @sort, vec(crtsec), 'descend'));
        elseif isinf(cutoffs(2, i))
            secs0{i} = crtsec(getnth(2, @sort, vec(crtsec), 'ascend'));
        else
            avg = mean(cutoffs(:, i));
            secs0{i} = crtsec(getnth(2, @sort, abs(vec(crtsec) - avg), 'descend'));
        end
    else
        if params.tails(i) > 1
            vec_mod = abs(vec - mu);
        else
            vec_mod = vec;
        end
        [~, order] = sort(vec_mod, 'descend');
        secs0{i} = order(1:params.cutoff(i));
        cutoffs(2, i) = vec_mod(secs0{i}(end));
        if params.tails(i) > 1
            cutoffs(1, i) = mu - cutoffs(2, i);
        end
    end

    % annotate the histograms
    if params.figures
        y_lim = get(gca, 'ylim');
        line([cutoffs(2, i) cutoffs(2, i)], y_lim, 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--');
        txt = ['>' num2str(cutoffs(2, i), '%0.2f')];
        text(cutoffs(2, i), 0.97*(y_lim(2) - y_lim(1)), txt, 'FontWeight', 'bold', 'FontSize', 11);
        if params.tails(i) > 1
            line([cutoffs(1, i) cutoffs(1, i)], y_lim, 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--');
            txt1 = ['<' num2str(cutoffs(1, i), '%0.2f')];
            
            text(cutoffs(1, i), 0.97*(y_lim(2) - y_lim(1)), txt1, 'FontWeight', 'bold', 'FontSize', 11);
        end
        xlabel(['comp #' num2str(i)], 'FontSize', 12, 'FontWeight', 'b');
        ylabel('Probability Density', 'FontSize', 12, 'FontWeight', 'b');
    end
end

% enforce non-overlapping condition, if required
if ~params.overlaps
    secs_orig = secs0;
    for i = 1:length(secs0)
        tmp = secs_orig([1:i-1 i+1:end]);
        other_secs = cell2mat(tmp(:));
        secs0{i} = setdiff(secs_orig{i}, other_secs);
    end
end

% output extra information, related to amount of 'sector'-ness
% this is currently non-public API
secdefs = cell(nfit, 1);
for i = 1:nfit
    if params.tails(i) == 1
        secdefs{i} = max(comps(:, i) - cutoffs(2, i), 0);
    else
        mu = (cutoffs(1, i) + cutoffs(2, i)) / 2;
        secdefs{i} = zeros(size(comps, 1), 1);
        secdefs{i}(secs0{i}) = abs(comps(secs0{i}, i) - mu);
    end
    % normalize
    secdefs{i} = secdefs{i}/max(secdefs{i});
end

% finally, make scatter plots of evecs vs. each other
if params.figures
    if nfit > 1
        nfigs = min(ceil(nfit/2), 4);
        
        % this allows the individual figures to have reasonable aspect ratios
        h = figure;
        set(h, 'Units', 'normalized', 'Position', [0 0.1 0.25*nfigs 0.3], 'Name', 'Sector Definition');
        
        % color the sectors according to the color wheel
        sector_hues = (0:nfit-1)/nfit;
        n = size(comps, 1);
        
        % default color: white
        res_hsv = repmat([1 0 1], n, 1);
        % keep track of residues that don't belong to any sector
        no_sector = 1:n;
        % when 'shade' is true, this value shows what fraction of the
        % points to be shaded; approximately 1 - fraction of the residues
        % will receive the full color of the sector they're part of
        fraction = 0.5;
        for i = 1:nfit
            if ~params.shade
                res_hsv(secs0{i}, 1) = sector_hues(i);
                res_hsv(secs0{i}, 2) = 1;
            else
                res_hsv(secs0{i}, 1) = sector_hues(i);
                res_hsv(secs0{i}, 2) = min(secdefs{i}(secs0{i}) / fraction, 1);
            end
            
            no_sector = setdiff(no_sector, secs0{i});
        end
        res_colors = hsv2rgb(res_hsv);
        
        for i = 1:2:nfit
            j = i + 1;
            if j > nfit
                j = 1;
            end
            subplot(1, nfigs, (i + 1)/2);

            smartscatter(comps(:, i), comps(:, j), 'shapes', '.', 'sizes', 25, ...
                'colors', res_colors);
            hold on;
            smartscatter(comps(no_sector, i), comps(no_sector, j), 'shapes', 'o', 'sizes', 5, ...
                'colors', 'k', 'tight', false);
            xlabel(['comp #' num2str(i)], 'FontSize', 12, 'FontWeight', 'b');
            ylabel(['comp #' num2str(j)], 'FontSize', 12, 'FontWeight', 'b');
            title(params.title, 'FontSize', 12);
            grid on;
        end
    end
end

extra = secdefs;
secs.secs = secs0;
secs.comps = comps;

% map to reference sequence
if ~isempty(params.refseq)
    alphawidths = arrayfun(@(s) length(s.map), params.refseq);
    mapped = cell(1, length(alphawidths));
    crt = 0;
    for i = 1:length(alphawidths)
        mapped{i} = cellfun(...
            @(s) params.refseq(i).map(...
                s(s > crt & s <= crt + alphawidths(i)) - crt), secs.secs, ...
                'uniform', false);
        crt = crt + alphawidths(i);
    end
    
    secs.secs_mapped = mapped;
end

end