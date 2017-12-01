function res = zeroblockdiag(m, structure, varargin)
% ZEROBLOCKDIAG Clear diagonal and near-diagonal blocks of a block matrix.
%   res = ZEROBLOCKDIAG(m, structure) clears the elements on the block
%   diagonal of the matrix m. The block structure is chosen to match the
%   'alphabets' and 'alphawidths' fields in 'structure'.
%
%   ZEROBLOCKDIAG(m, structure, d) removes elements up to the dth diagonal
%   (i.e., d = 0 is only the main diagonal; d = 1 removes the elements one
%   above and one below the main diagonal; etc.).
%
%   ZEROBLOCKDIAG(m, L, ...) assumes that all the blocks are L x L. This is
%   equivalent to having a structure with a single alphabet containing L
%   letters (except the gap).
%
%   Options:
%    'value' <x>
%       The value with which to fill the cleared parts of the matrix.
%       (default: 0)

% Tiberiu Tesileanu (2012-2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('d', 0, @(x) isnumeric(x) && isscalar(x));

parser.addParamValue('value', 0, @(x) isnumeric(x) && isscalar(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

if size(m, 1) ~= size(m, 2)
    error([mfilename ':badsize'], 'The matrix should be square.');
end

res = m;

if ~isstruct(structure)
    if ~isnumeric(structure) || ~isscalar(structure)
        error([mfilename ':badstruc'], 'The second argument should either be an alignment-like structure or an integer.');
    end
    L = structure;
    if mod(size(m, 1), L) ~= 0
        error([mfilename ':sizemmatch'], 'Mismatch between size of m and structure.');
    end
    Npos = size(m, 1) / L;
    for i = 1:Npos
        iidxs = L*(i-1) + (1:L);
        res(iidxs, iidxs) = params.value;
        for j = (i+1):min(i+params.d, Npos)
            jidxs = L*(j-1) + (1:L);
            
            res(iidxs, jidxs) = params.value;
            res(jidxs, iidxs) = params.value;
        end
    end
else
    binmap = getbinmap(structure);
    expectedsize = binmap{end}(end);
    if expectedsize ~= size(m, 1) || expectedsize ~= size(m, 2)
        error([mfilename ':sizemmatch'], 'Mismatch between size of m and structure.');
    end

    Npos = length(binmap);
    
    for i = 1:Npos
        iidxs = binmap{i};
        res(iidxs, iidxs) = params.value;
        for j = (i + 1):min(i + params.d, Npos)
            jidxs = binmap{j};
            
            res(iidxs, jidxs) = params.value;
            res(jidxs, iidxs) = params.value;
        end
    end
end

end