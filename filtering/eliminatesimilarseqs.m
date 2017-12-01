function newalign = eliminatesimilarseqs(alignment, max_seqid, varargin)
% ELIMINATESIMILARSEQS Eliminate repeated or nearly-repeated sequences.
%   align = ELIMINATESIMILARSEQS(alignment, max_seqid) eliminates sequence
%   repetitions, where repetition is defined as sequence identity greater
%   than or equal to max_seqid. This sets the sequence weights to
%   identically 1.
%
%   In more detail, the function creates an undirected graph with sequences
%   as vertices, in which an edge between two sequences means that they are
%   identical within the max_seqid threshold. The resulting alignment is
%   created by choosing a representative for each of the connected
%   components of this graph. The representative is chosen to have the
%   samllest number of gaps; if there are several such sequences, an
%   arbitrary choice is made.

% Tiberiu Tesileanu (2012-2014)

% allow other functions to check whether this function uses a C++ module
if nargin == 1 && ischar(alignment) && strcmp(alignment, 'hascpp')
    newalign = (exist('calcseqw_cpp', 'file') == 3);
    return;
end

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('nocpp', false, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

if ~alncheck(alignment)
    error([mfilename ':badarg'], 'The first argument should an alignment structure.');
end

% nothing to do if max_seqid is > 1
if max_seqid > 1
    newalign = alignment;
    return;
end

% the mask is the adjacency matrix for the graph
if ~params.nocpp && exist('calcseqw_cpp', 'file') == 3
    mask = calcseqw_cpp(alignment.data, max_seqid);
else
    numdata = double(alignment.data);
    mask = sparse(squareform(pdist(numdata, 'hamm') <= 1 - max_seqid));
end

% find the connected components
[n_comps, comp_labels] = graphconncomp(mask, 'Directed', false);
    
% now we know what size the resulting alignment will be
newdata = repmat(blanks(size(alignment.data, 2)), n_comps, 1);
newannot = cell(n_comps, 1);

% get gaps per sequence
gapstruct = getgapstructure(alignment);
gapsperseq = sum(gapstruct, 2);

% go through the components and find the representatives
for i = 1:n_comps
    selection = (comp_labels == i);
    seqs = alignment.data(selection, :);
    
    % choose the sequence with the fewest gaps (the first one if there are
    % several)
    % (if there is a focal sequence, this should choose it automatically,
    % since it had no gaps, and it was first in the alignment)
    [~, indices] = min(gapsperseq(selection));
    newdata(i, :) = seqs(indices(1), :);
    annots = alignment.annotations(selection);
    newannot{i} = annots{indices(1)};
end

% make sure to copy all other members
newalign = alignment;

% update the fields
newalign.data = newdata;
newalign.seqw = ones(n_comps, 1);
newalign.annotations = newannot;

end