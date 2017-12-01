function [hvals_all, hvals_sel] = makesinglehist(n, comps, values_exp, compmap, expmap, varargin)
% MAKESINGLEHIST Make a histogram comparing experimental data within a
% group to the entire data set.
%   [hvals_all, hvals_sel] = MAKESINGLEHIST(n, comps, values_exp, compmap, expmap)
%   makes a histogram of the experimental values given in values_exp, and
%   overlays on top of it a histogram with the values corresponding to the
%   residues that have the n largest values in comps. The mapping between
%   the components of comps and the experimental measurements is done using
%   the maps expmap and compmap. Expmap(i) gives an index in some reference
%   sequence corresponding to values_exp(i); compmap(i) gives the index in
%   the same reference sequence that corresponds to comps(i).
%
%   The function returns the two sets of values that were histogrammed.
%
%   Options:
%    'figures' <b>
%       Whether to draw the actual histogram. Set to false to only return
%       hvals_all and hvals_sel, without any histogram being displayed.
%       (default: true)
%    'range' <x>
%       Bin centers to use for the histograms, or number of bins.
%       (default: 25)
%    'title' <s>
%       String to use for the figure's title.
%       (default: none)

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('figures', true, @(b) islogical(b) && isscalar(b));
parser.addParamValue('range', 25, @(v) isnumeric(v) && isvector(v));
parser.addParamValue('title', '', @(s) ischar(s) && isvector(s));

% parse the options
parser.parse(varargin{:});
params = parser.Results;

[idxs0_comps, idxs0_exp] = findcommonidxs(compmap, expmap);
hvals_all = values_exp(idxs0_exp);
[~, sorted_comps] = sort(comps, 'descend');
sel_comps = sorted_comps(1:n);
hvals_sel = hvals_all(ismember(idxs0_comps, sel_comps));

if params.figures
    multihist({hvals_all, hvals_sel}, params.range, ...
        'mode', 'overlap', ...
        'colors', [[0.7 0.7 0.7] ; [0 0 0]], ...
        'width', 0.9);
    tmp = axis;
    axis([min(params.range) + 0.05 max(params.range) + 0.05 tmp(3) tmp(4)]);
    ylabel('counts');
    xlabel('mutational effect <\Delta E_i^x>_x');
    if ~isempty(params.title)
        title(params.title);
    end
    beautifygraph;
    set(gca, 'xgrid', 'on', 'ygrid', 'on', ...
        'xminorgrid', 'off', 'yminorgrid', 'off');
end

end