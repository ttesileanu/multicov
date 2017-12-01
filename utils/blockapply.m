function res = blockapply(mat, fct, pattern, varargin)
% BLOCKAPPLY Apply a function to each block of a mixed-alphabet block
% matrix or vector.
%   res = BLOCKAPPLY(mat, fct, pattern) applies the function fct to each
%   block of the vector or matrix mat. The block structure is inferred from
%   the alignment-like structure 'pattern' (which should have fields
%   'alphabets' and 'alphawidths').
%
%   The result is typically placed in a vector or a matrix. For this to
%   work, each application of fct must return a scalar. Use the
%   'symmetric', true flag to optimize for cases in which the input is a
%   matrix, and the result is expected to be symmetric.
%
%   If the first application of fct returns a non-scalar, the output is
%   placed in a cell array instead. This behavior can be forced if
%   'uniform', true is used.
%
%   Options:
%    'indices' <b>
%       If true, the function 'fct' is also passed the index/indices of the
%       block it's acting on (i.e., fct(block, i) for vector input, and
%       fct(block, i, j) for matrix input).
%       (default: false)
%    'symmetric' <b>
%       Whether to assume the result is symmetric. This only applies to
%       matrix inputs. When the output is a cell array, res(j, i) is set to
%       res(i, j)' (i.e., it is assumed that a transpose is needed).
%       (default: true)
%    'uniform' <b>
%       Whether to assume each application of fct will return a scalar.
%       (default: autodetect based on the first function call for fct)

% Tiberiu Tesileanu (2012, 2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('indices', false, @(b) isscalar(b) && islogical(b));
parser.addParamValue('symmetric', true, @(b) isscalar(b) && islogical(b));
parser.addParamValue('uniform', [], @(b) isscalar(b) && islogical(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

binmap = getbinmap(pattern);
% autodetect uniformity, if required
if isempty(params.uniform)
    % calling fct one more time than needed, but seems worth it rather than
    % having ifs in the loops
    if isvector(mat)
        if params.indices
            tmp = fct(mat(binmap{1}), 1);
        else
            tmp = fct(mat(binmap{1}));
        end
    else
        if params.indices
            tmp = fct(mat(binmap{1}, binmap{1}), 1, 1);
        else
            tmp = fct(mat(binmap{1}, binmap{1}));
        end
    end
    
    params.uniform = isscalar(tmp);
end

n = length(binmap);
ressz = [n, n];
% conserve shape of input
ressz(size(mat) == 1) = 1;

if params.uniform
    res = zeros(ressz);
else
    res = cell(ressz);
end
if n == 0
    return;
end

% handling each case separately for efficiency and clarity
% an unfortunate consequence is that there is *a lot* of code duplication
if params.indices
    if isvector(mat)
        if params.uniform
            for i = 1:n
                res(i) = fct(mat(binmap{i}), i);
            end
        else
            for i = 1:n
                res{i} = fct(mat(binmap{i}), i);
            end
        end
    else
        if params.uniform
            if params.symmetric
                for i = 1:n
                    idxsi = binmap{i};
                    for j = (i+1):n
                        idxsj = binmap{j};
                        res(i, j) = fct(mat(idxsi, idxsj), i, j);
                    end
                end
                res = res + res';
                for i = 1:n
                    idxsi = binmap{i};
                    res(i, i) = fct(mat(idxsi, idxsi), i, i);
                end
            else
                for i = 1:n
                    idxsi = binmap{i};
                    for j = 1:n
                        idxsj = binmap{j};
                        res(i, j) = fct(mat(idxsi, idxsj), i, j);
                    end
                end
            end
        else
            if params.symmetric
                for i = 1:n
                    idxsi = binmap{i};
                    res{i, i} = fct(mat(idxsi, idxsi), i, i);
                    for j = (i+1):n
                        idxsj = binmap{j};
                        res{i, j} = fct(mat(idxsi, idxsj), i, j);
                        res{j, i} = res{i, j}';
                    end
                end
            else
                for i = 1:n
                    idxsi = binmap{i};
                    for j = 1:n
                        idxsj = binmap{j};
                        res{i, j} = fct(mat(idxsi, idxsj), i, j);
                    end
                end
            end
        end
    end
else
    if isvector(mat)
        if params.uniform
            for i = 1:n
                res(i) = fct(mat(binmap{i}));
            end
        else
            for i = 1:n
                res{i} = fct(mat(binmap{i}));
            end
        end
    else
        if params.uniform
            if params.symmetric
                for i = 1:n
                    idxsi = binmap{i};
                    for j = (i+1):n
                        idxsj = binmap{j};
                        res(i, j) = fct(mat(idxsi, idxsj));
                    end
                end
                res = res + res';
                for i = 1:n
                    idxsi = binmap{i};
                    res(i, i) = fct(mat(idxsi, idxsi));
                end
            else
                for i = 1:n
                    idxsi = binmap{i};
                    for j = 1:n
                        idxsj = binmap{j};
                        res(i, j) = fct(mat(idxsi, idxsj));
                    end
                end
            end
        else
            if params.symmetric
                for i = 1:n
                    idxsi = binmap{i};
                    res{i, i} = fct(mat(idxsi, idxsi));
                    for j = (i+1):n
                        idxsj = binmap{j};
                        res{i, j} = fct(mat(idxsi, idxsj));
                        res{j, i} = res{i, j}';
                    end
                end
            else
                for i = 1:n
                    idxsi = binmap{i};
                    for j = 1:n
                        idxsj = binmap{j};
                        res{i, j} = fct(mat(idxsi, idxsj));
                    end
                end
            end
        end
    end
end

end