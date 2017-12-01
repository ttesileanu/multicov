function [newalign, indices] = filtercolumns(alignment, max_gaps, varargin)
% FILTERCOLUMNS Filter alignment positions that have too many gaps.
%   newalign = FILTERCOLUMNS(alignment, max_gaps) removes those
%   positions in the alignment that have a fraction of gaps larger than
%   max_gaps.
%
%   FILTERCOLUMNS(..., 'edges', true) only removes columns from the
%   edges of the alignment. That is, it removes the left-most columns with
%   gap fraction larger than max_gap, and the rigt-most columns with the
%   same property. Columns in the middle that might also have many gaps are
%   kept.

% Tiberiu Tesileanu (2012-2014)

if ~alncheck(alignment)
    error([mfilename ':badarg'], 'The first argument should an alignment structure.');
end

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('edges', false, @(b) isscalar(b) && islogical(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

% find the gaps
gaps = getgapstructure(alignment);

% get rid of positions with too many gaps
Neff = sum(alignment.seqw);
gapfraction = alignment.seqw(:)'*gaps/Neff;
mask1 = (gapfraction <= max_gaps);

if params.edges
    % need to update the mask so that it only has zeros on the edges
    lidx = find(mask1, 1, 'first');
    ridx = find(mask1, 1, 'last');
    if ~isempty(lidx)
        mask1(lidx:ridx) = true;
    end
end

newdata = alignment.data(:, mask1);

% more involved: we need to update the alphabets
newalphabets = {};
newalphawidths = [];
ranges = getalpharanges(alignment);
indices = cell(length(alignment.alphabets), 1);
for i = 1:length(alignment.alphabets)
    submask = mask1(ranges(1, i):ranges(2, i));
    indices{i} = find(submask);
    indices{i} = indices{i}(:);
    newwidth = length(indices{i});
    if newwidth > 0
        newalphabets = [newalphabets alignment.alphabets{i}]; %#ok<AGROW>
        newalphawidths = [newalphawidths newwidth]; %#ok<AGROW>
    end
end

% make sure to copy all the fields from the input structure
newalign = alignment;

newalign.alphabets = newalphabets;
newalign.alphawidths = newalphawidths;
newalign.data = newdata;

% update mapping to reference sequences
crtj = 1;
for i = 1:length(alignment.alphabets)
    if ~isempty(indices{i})
        newalign.refseq(crtj).map = alignment.refseq(i).map(indices{i});
        crtj = crtj + 1;
    end
end
nonemptymask = ~cellfun(@isempty, indices);
newalign.refseq = newalign.refseq(nonemptymask);

end