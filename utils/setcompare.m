function sim = setcompare(set1, set2, varargin)
% SETCOMPARE Compare two sets.
%   sim = SETCOMPARE(set1, set2) returns a measure of the similarity of two
%   sets. By default this is defined as the ratio between the number of
%   elements in the intersection of set1 and set2, and the mean of the
%   numbers of elements in the sets.
%
%   Note that set1 and set2 may either be arrays of numbers or cell arrays
%   of strings. If they are of any other type, or if the type of set1 is
%   different from that of set2, the functionally immediately returns 0.
%
%   Options:
%    'method' <s>
%       Choose the method for calculating the similarity. This is always
%       the ratio between the number of elements in the intersection of
%       set1 and set2, and one of the following:
%        'union':   number of elements in union of set1 and set2
%        'min':     number of elements in smaller set
%        'max':     number of elements in lager set
%        'mean':    average number of elements between set1 and set2
%        'first':   number of elements in set1
%        'second':  number of elements in set2
%       (default: 'mean')

% Tiberiu Tesileanu (2012, 2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('method', 'mean', @(s) ismember(s, {'union', 'min', 'max', 'mean', 'first', 'second'}));

% parse
parser.parse(varargin{:});
params = parser.Results;

set1 = set1(:);
set2 = set2(:);
if iscell(set1)
    if ~iscell(set2)
        sim = 0;
        return;
    elseif (~isempty(set1) && ~ischar(set1{1})) || (~isempty(set2) && ~ischar(set2{1}))
        sim = 0;
        return;
    end
elseif iscell(set2)
    sim = 0;
    return;
elseif ~isnumeric(set1) || ~isnumeric(set2)
    sim = 0;
    return;
end;

nint = numel(intersect(set1, set2));

switch params.method
    case 'union'
        ntot = numel(union(set1, set2));
    case 'min'
        ntot = min(numel(set1), numel(set2));
    case 'max'
        ntot = max(numel(set1), numel(set2));
    case 'mean'
        ntot = (numel(set1) + numel(set2)) / 2;
    case 'first'
        ntot = numel(set1);
    case 'second'
        ntot = numel(set2);
end

sim = nint/ntot;

end