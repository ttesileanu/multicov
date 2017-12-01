function dist = getdistmat(coords, varargin)
% GETDISTMAT Generate a distance matrix given atom positions or PDB
% structure.
%   dist = GETDISTMAT(coords) returns a matrix dist such that dist(i, j) is
%   the distance between residues i and j, given a cell array in which each
%   cell gives the atom positions for one residue, represented as a
%   structure with fields 'pos', a K x 3 matrix in which each row is the
%   position of one atom, and 'type', a cell array of K strings giving the
%   types of the atoms (this is the structure returned by getpdbcoords).
%   By default, the distance returned is the minimal inter-atom distance.
%
%   dist = GETDISTMAT({pdb, chain}) uses the coordinates from the given
%   chain of the given PDB structure.
%
%   Options:
%    'method' <s/f>
%       Either a string identifying the kind of distance calculated for
%       each pair of residues,
%         'min':    minimal inter-atom distance
%         'gcm':    center of mass distance, assuming all atoms are the
%                   same weight
%         'ca':     distance between CA atoms
%       or a function handle taking the two structures for the two
%       residues, and returning the distance.
%       (default: 'min')
%
% See also: GETPDBCOORDS.

% Tiberiu Tesileanu (2012-2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('method', 'min', @(x) ischar(x) || isa(x, 'function_handle'));

% parse
parser.parse(varargin{:});
params = parser.Results;

if iscell(coords) && length(coords) == 2 && isstruct(coords{1}) && ischar(coords{2})
    coords = getpdbcoords(coords{1}, coords{2});
end

switch params.method
    case 'min'
        params.method = @mindist;
    case 'gcm'
        params.method = @gcmdist;
    case 'ca'
        params.method = @cadist;
end

n = length(coords);
dist = zeros(n);
for i = 1:n
    for j = (i + 1):n
        dist(i, j) = params.method(coords{i}, coords{j});
    end
end
dist = dist + dist';

end

function d = mindist(r1, r2)

distances = pdist2(r1.pos, r2.pos);
d = min(distances(:));

end

function d = gcmdist(r1, r2)

d = norm(mean(r1.pos) - mean(r2.pos));

end

function d = cadist(r1, r2)

mask1 = find(strcmp(r1.type, 'CA'), 1);
mask2 = find(strcmp(r2.type, 'CA'), 1);
if isempty(mask1) || isempty(mask2)
    error([mfilename ':noca'], 'Missing CA atom.');
end
d = norm(r1.pos(mask1, :) - r2.pos(mask2, :));

end