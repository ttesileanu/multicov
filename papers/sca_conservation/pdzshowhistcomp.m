function pdzshowhistcomp(expresults, alnmap, seccomps, cons, scadiag, varargin)
% PDZSHOWHISTCOMP A helper for the pdzanalysis script that creates
% histograms similar to Figure 3 from McLaughlin Jr. et al. 2012.
%   PDZSHOWHISTCOMP(expresults, pdbats, seccomps, cons, scadiag)
%   makes a figure with three histograms, showing how different the effect
%   of mutating sector positions, top conservation positions, or top
%   diag(SCA) positions is to mutating arbitrary positions. By default,
%   each of the histograms contains 20 positions. The mapping between the
%   experimental results and the alignment data is done using alnmap and
%   expresults.pos, and the experimental data used is expresults.Eavgx.
%
%   If any of seccomps, cons, or scadiag are empty, the corresponding
%   histogram is not drawn.
%
%   Options:
%    'comment' <s>
%       Text to add to the titles of figures.
%       (default: none)
%    'displayks' <b>
%       Whether to display the Kolmogorov-Smirnov p-values for comparisons
%       between the different histograms.
%       (default: true)
%    'figures' <b>
%       Whether to draw the histograms. Set to false to only display the
%       p-values for the statistical tests.
%       (default: true)
%    'format' <n/s>
%       Number format to use for displaying the p-values (sprintf format),
%       or an integer showing the number of significant digits (as for
%       num2str).
%       (default: 2)
%    'n' <n>
%       The number of top positions to consider from seccomps, cons, and
%       scadiag.
%       (default: 20)
%    'pvals' <b>
%       Whether to show the p-values on the histograms.
%       (default: false)
%    'titles' {<s1>, <s2>, <s3>}
%       Titles to use for the plots.
%       (default: {'Sector positions', 'Top conserved positions',
%        'Top positions by diag(SCA)'})
%
% See also: MAKESINGLEHIST.

% Tiberiu Tesileanu (2013-2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('comment', '', @(s) ischar(s) && isvector(s));
parser.addParamValue('displayks', true, @(b) islogical(b) && isscalar(b));
parser.addParamValue('figures', true, @(b) islogical(b) && isscalar(b));
parser.addParamValue('format', 2, @(x) (isscalar(x) && isnumeric(x)) || (ischar(x) && isvector(x)));
parser.addParamValue('n', 20, @(x) isnumeric(x) && isscalar(x));
parser.addParamValue('pvals', false, @(b) isscalar(b) && islogical(b));
parser.addParamValue('titles', ...
    {'Sector positions', ...
     'Top conserved positions', ...
      'Top positions by diag(SCA)'}, ...
    @(c) iscell(c) && isvector(c) && length(c) == 3 && all(cellfun(@(s) isvector(s) && ischar(s), c)));

% parse options

parser.parse(varargin{:});
params = parser.Results;

%[mask0_aln, mask0_exp] = findcommonmasks(alnpdbats, expresults.pos_pdb);
%values0_exp = expresults.Eavgx(mask0_exp);

if params.figures
    figure;
end

plotmask = cellfun(@(c) ~isempty(c), {seccomps, cons, scadiag});
nplots = sum(plotmask);
crt = 1;

if plotmask(1)
    % how the sector correlates with functional significance
    if params.figures
        subplot(nplots, 1, crt);
    end
    crt = crt + 1;
    crttitle = params.titles{1};
    if ~isempty(params.comment)
        crttitle = [crttitle ' ' params.comment];
    end
    crttitle = [crttitle ' (' int2str(params.n) ' residues)'];
    [hvals_sec_all, hvals_sec_sel] = makesinglehist(params.n, seccomps, expresults.Eavgx, ...
        alnmap, expresults.pos, 'title', crttitle, 'range', (-1.5:0.1:0.5) - 0.05, ...
        'figures', params.figures);
end

if plotmask(2)
    % how conservation correlates with functional significance
    if params.figures
        subplot(nplots, 1, crt);
    end
    crt = crt + 1;
    crttitle = params.titles{2};
    if ~isempty(params.comment)
        crttitle = [crttitle ' ' params.comment];
    end
    crttitle = [crttitle ' (' int2str(params.n) ' residues)'];
    [hvals_cons_all, hvals_cons_sel] = makesinglehist(params.n, cons, expresults.Eavgx, ...
        alnmap, expresults.pos, 'title', crttitle, 'range', (-1.5:0.1:0.5) - 0.05, ...
        'figures', params.figures);
end

if plotmask(3)
    % how diag(SCA) correlates with functional significance
    if params.figures
        subplot(nplots, 1, crt);
    end
    crttitle = params.titles{3};
    if ~isempty(params.comment)
        crttitle = [crttitle ' ' params.comment];
    end
    crttitle = [crttitle ' (' int2str(params.n) ' residues)'];
    [hvals_diag_all, hvals_diag_sel] = makesinglehist(params.n, scadiag, expresults.Eavgx, ...
        alnmap, expresults.pos, 'title', crttitle, 'range', (-1.5:0.1:0.5) - 0.05, ...
        'figures', params.figures);
end

if params.figures
    crttitle = 'Efficacy of SCA over single-site methods';
    if ~isempty(params.comment)
        crttitle = [crttitle ' ' params.comment];
    end
    set(gcf, 'units', 'normalized', 'position', [0.35 0.1 0.3 0.8*nplots/3], ...
        'name', crttitle);
end

% how different are the distributions?
if params.displayks
%    disp('Statistical two-sample tests:');
    
    crt = 1;
    textopts = {'FontName', 'Helvetica', 'FontSize', 12};
    if plotmask(1)
        crttitle = [params.titles{1} ' vs. all positions (n = ' int2str(params.n) ')'];
        if ~isempty(params.comment)
            crttitle = [crttitle ' ' params.comment];
        end
        disp(crttitle);
        [~, ksp] = kstest2(hvals_sec_sel, hvals_sec_all, 0.05, 'larger');
        mwwp = ranksum(hvals_sec_sel, hvals_sec_all, 'tail', 'left');
        disp(['Kolmogorov-Smirnov p = ' num2str(ksp, params.format) ', Mann-Whitney p = ' num2str(mwwp, params.format)]);
        
        if params.pvals && params.figures
            subplot(nplots, 1, crt);
            tmp = axis;
            text(-1.45, 0.25*tmp(3) + 0.75*tmp(4), ['Mann-Whitney p = ' num2str(mwwp, params.format)], ...
                textopts{:});
            crt = crt + 1;
        end
    end

    if plotmask(2)
        crttitle = [params.titles{2} ' vs. all positions (n = ' int2str(params.n) ')'];
        if ~isempty(params.comment)
            crttitle = [crttitle ' ' params.comment];
        end
        disp(crttitle);
        [~, ksp] = kstest2(hvals_cons_sel, hvals_cons_all, 0.05, 'larger');
        mwwp = ranksum(hvals_cons_sel, hvals_cons_all, 'tail', 'left');
        disp(['Kolmogorov-Smirnov p = ' num2str(ksp, params.format) ', Mann-Whitney p = ' num2str(mwwp, params.format)]);
        
        if params.pvals && params.figures
            subplot(nplots, 1, crt);
            tmp = axis;
            text(-1.45, 0.25*tmp(3) + 0.75*tmp(4), ['Mann-Whitney p = ' num2str(mwwp, params.format)], ...
                textopts{:});
            crt = crt + 1;
        end
        
        crttitle = [params.titles{1} ' vs. ' params.titles{2} ' (n = ' int2str(params.n) ')'];
        if ~isempty(params.comment)
            crttitle = [crttitle ' ' params.comment];
        end
        disp(crttitle);
        [~, ksp] = kstest2(hvals_sec_sel, hvals_cons_sel, 0.05);
        mwwp = ranksum(hvals_sec_sel, hvals_cons_sel);
        disp(['Kolmogorov-Smirnov p = ' num2str(ksp, params.format) ', Mann-Whitney p = ' num2str(mwwp, params.format)]);
    end

    if plotmask(3)
        crttitle = [params.titles{3} ' vs. all positions (n = ' int2str(params.n) ')'];
        if ~isempty(params.comment)
            crttitle = [crttitle ' ' params.comment];
        end
        disp(crttitle);
        [~, ksp] = kstest2(hvals_diag_sel, hvals_diag_all, 0.05, 'larger');
        mwwp = ranksum(hvals_diag_sel, hvals_diag_all, 'tail', 'left');
        disp(['Kolmogorov-Smirnov p = ' num2str(ksp, params.format) ', Mann-Whitney p = ' num2str(mwwp, params.format)]);
        
        if params.pvals && params.figures
            subplot(nplots, 1, crt);
            tmp = axis;
            text(-1.45, 0.25*tmp(3) + 0.75*tmp(4), ['Mann-Whitney p = ' num2str(mwwp, params.format)], ...
                textopts{:});
        end
     
        crttitle = [params.titles{1} ' vs. ' params.titles{3} ' (n = ' int2str(params.n) ')'];
        if ~isempty(params.comment)
            crttitle = [crttitle ' ' params.comment];
        end
        disp(crttitle);
        [~, ksp] = kstest2(hvals_sec_sel, hvals_diag_sel, 0.05);
        mwwp = ranksum(hvals_sec_sel, hvals_diag_sel);
        disp(['Kolmogorov-Smirnov p = ' num2str(ksp, params.format) ', Mann-Whitney p = ' num2str(mwwp, params.format)]);
    end
end

if params.figures
    preparegraph;
end

end