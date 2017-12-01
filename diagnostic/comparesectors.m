function [sim, map, mat] = comparesectors(secs1, secs2, varargin)
% COMPARESECTORS Compare groups of sectors.
%   [sim, map, mat] = COMPARESECTORS(secs1, secs2) compares the two sets of
%   sectors. Each of secs1 and secs2 should be a cell array of arrays
%   (either numeric or cell), each of which identifies a sector. If either
%   set contains a single sector, it can be provided without the outer cell
%   array. Note that for reference-sequence-mapped sectors, COMPARESECTORS
%   has to be run separately for each alphabet.
%
%   The sectors are compared as sets of elements, using the setcompare
%   function. The order of the sectors is allowed to be different; the
%   optimal mapping between sectors in secs1 and sectors in secs2 is found
%   using a version of the Hungarian algorithm. The implementation,
%   assignmentoptimal.m, was downloaded from
%   http://www.mathworks.com/matlabcentral/fileexchange/6543
%
%   The optimal mapping is returned in 'map': secs1{i} is most similar to
%   secs2{map(i)}; the level of similarity is returned in sim(i); and the
%   matrix 'mat' has elements mat(i, j) = similarity between secs1(i) and
%   secs2(j).
%
%   Options:
%    'method' <s>
%       How to caculate set similarity. See setcompare for the possible
%       options.
%       (default: same as setcompare)
%
% See also SETCOMPARE.

% Tiberiu Tesileanu (2012, 2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('method', '', @(s) ischar(s) && isvector(s));

% parse
parser.parse(varargin{:});
params = parser.Results;

% need to handle the various input formats
if ~iscell(secs1) || isempty(secs1) || ischar(secs1{1})
    secs1 = {secs1};
end
if ~iscell(secs2) || isempty(secs2) || ischar(secs2{1})
    secs2 = {secs2};
end

nsecs1 = length(secs1);
nsecs2 = length(secs2);

% calculate the matrix of pairwise similarities
if isempty(params.method)
    methargs = {};
else
    methargs = {'method', params.method};
end
sim_matrix = zeros(nsecs1, nsecs2);
for i = 1:nsecs1
    for j = 1:nsecs2
        sim_matrix(i, j) = setcompare(secs1{i}, secs2{j}, methargs{:});
    end
end

% find the optimal matching
sim = zeros(nsecs1, 1);
map = assignmentoptimal(max(sim_matrix(:)) - sim_matrix);
for i = 1:nsecs1
    if map(i) > 0
        sim(i) = sim_matrix(i, map(i));
    end
end
mat = sim_matrix;

end