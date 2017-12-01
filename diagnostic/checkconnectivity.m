function [conn, info] = checkconnectivity(positions, contacts, varargin)
% CHECKCONNECTIVITY Check to what extent a set of positions forms a
% connected structural component.
%   [conn, info] = CHECKCONNECTIVITY(positions, contacts) checks whether
%   the set of positions forms a connected structural component. The matrix
%   'contacts' is an adjacency matrix (e.g., a contact map for a tertiary
%   structure of some molecule), and 'positions' is a vector of indices.
%   The measure of connectivity employed is the ratio between the size of
%   the largest connected component, and the total size of the 'positions',
%   rescaled and shifted such that it is 0 when the largest component is
%   made up of a single node, and 1 when the largest component comprises
%   all of the positions. The structre 'info' contains a field 'compsizes'
%   that returns the sizes of all the connected components, in descending
%   order.
%
%   CHECKCONNECTIVITY(..., 'pvalue', true) also returns an estimate of the
%   p-value for the null hypothesis that the observed connectivity of the
%   set of positions could have been obtained by choosing a random set of
%   positions with the same size. This is returned as a field 'p' of the
%   'info' structure. Note that the value is estaimted from many random
%   samples, and as such it can change between calls to this function.
%
%   Options:
%    'nsamples'
%       Number of samples to use when calculating the p-value.
%       (default: 1000)
%    'pvalue' <b>
%       If true, calculate and return the p-value.
%       (default: false)

% Tiberiu Tesileanu (2012, 2014)

% XXX something smart should be done about multi-alphabet cases...
% XXX I don't think the thing below is the best
% if iscell(contacts)
%     % dealing with multiple alphabets
%     split = alnsplitindices(cellfun(@(map) size(map, 1), contacts), positions);
%     n = length(contacts);
%     connres = cell(1, n);
%     compsizes = cell(1, n);
%     for i = 1:n
%         [connres{i}, compsizes{i}] = checkconnectivity(split{i}, contacts{i}, varargin{:});
%     end
%     
%     return;
% end

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('pvalue', false, @(b) isscalar(b) && islogical(b));
parser.addParamValue('nsamples', 1000, @(x) isscalar(x) && isnumeric(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

% handle the trivial cases
if isempty(positions) || length(positions) == 1
    conn = 1;
    info.compsizes = ones(length(positions), 1);
    if params.pvalue
        info.p = 1;
    end
    return;
end

adj = sparse(contacts(positions, positions));

[ncomps, assignments] = graphconncomp(adj, 'directed', false);

compsizes = zeros(1, ncomps);
for i = 1:ncomps
    compsizes(i) = sum(assignments == i);
end

compsizes = sort(compsizes, 'descend');
% shift and rescale, so that if the largest component has one element,
% connectivity = 0, if the whole group is connected, connectivity = 1
conn = (compsizes(1) - 1) / (length(positions) - 1);

info.compsizes = compsizes;

% calculate the p-value, if required
if params.pvalue
    rndconns = zeros(params.nsamples, 1);
    npos = length(positions);
    for i = 1:params.nsamples
        % get a random selection of residues
        positions = randperm(size(contacts, 1), npos);
        rndconns(i) = checkconnectivity(positions, contacts, 'pvalue', false);
    end
    info.p = sum(rndconns >= conn) / params.nsamples;
end

end