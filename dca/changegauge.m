function newparams = changegauge(params, gauge)
% CHANGEGAUGE Switch to an arbitrary gauge in a maximum entropy model.
%   newparams = CHANGEGAUGE(params, gauge) switches the maximum entropy
%   model couplings represented by the given parameters structure to an
%   arbitrary gauge. 'Gauge' is a cell array which gives, for each site i
%   in the alignment, a vector v_i. The gauge is defined by the constraints
%       J_{ij}*v_j = 0,   dot(h_i, v_i) = 0,    for all i, j,
%   where J_{ij} is the block in the coupling matrix corresponding to sites
%   i and j, and h_i is the block in the fields corresponding to site i.
%
%   Alternatively, the gauge can be provided as a single long vector, in
%   which all the cells in the cell array described above are concatenated.
%
%   The gauge can also be given for each of the alphabets in params, in
%   which case the number of elements in 'gauge' must be equal to this
%   number of alphabets, and each position corresponding to a give alphabet
%   uses the same gauge. For the single-alphabet case, 'gauge' can be the
%   gauge vector v itself instead of a cell array with this vector as the
%   single element.
%
%   newparams = CHANGEGAUGE(params, 'zerosum') is a special case of the
%   above, in which the gauge is the 'zero-sum' gauge, corresponding to
%   v_i(a) = 1 for all characters a, and for any i.
%
%   newparams = CHANGEGAUGE(params, 'gap') is a different special case in
%   which the gap entries of the couplings and fields are identically zero.
%   In other words, for all sites, v_i(gap) = 1, and v_i(a ~= gap) = 0.

% Tiberiu Tesileanu (2013-2014)

% need to identify which of the zillion forms of the function we are dealing
% with
binmap = getbinmap(params);
len = binmap{end}(end);
n = length(binmap);
if isnumeric(gauge)
    if length(params.alphabets) == 1 && length(gauge) == length(binmap{1})
        gauge = repmat(gauge(:), params.alphawidths, 1);
    end
    if length(gauge) ~= len
        error([mfilename ':badvlen'], 'The gauge vector input has invalid length.');
    end
    gauge = cellfun(@(idxs) gauge(idxs), binmap, 'uniform', false);
elseif iscell(gauge)
    if length(gauge) == length(params.alphabets)
        newgauge = {};
        for i = 1:length(gauge)
            newgauge = [newgauge ; repmat(gauge(i), params.alphawidths(i), 1)]; %#ok<AGROW>
        end
        gauge = newgauge;
    end
    if length(gauge) ~= n
        error([mfilename ':badclen'], 'The gauge cell array input has invalid length.');
    end
elseif ischar(gauge)
    switch gauge
        case 'zerosum'
            gauge = cellfun(@(idxs) ones(size(idxs)), binmap, 'uniform', false);
        case 'gap'
            gauge = cellfun(@(idxs) [1 ; zeros(length(idxs)-1, 1)], binmap, 'uniform', false);
        otherwise
            error([mfilename ':badgaugestr'], ['Unrecognized gauge choice, ''' gauge '''.']);
    end
else
    error([mfilename ':badgauge'], 'The gauge should be a numeric vector, a cell array, or a string.');
end

% by now 'gauge' should be a cell array with n elements, each of which is a
% numeric vector identifying the gauge

newparams = params;
newcouplings = zeros(size(params.couplings));
for i = 1:n
    idxsi = binmap{i};
    vi = gauge{i}(:);
    
    % this normalization makes the formulas simpler
    vi = vi/sum(vi);
    
    % the diagonal needs special treatment
    
    % the blocks on the diagonal are themselves diagonal, so work with only
    % the non-zero elements
    dblock0 = diag(params.couplings(idxsi, idxsi));
    dblock = dblock0 - dot(dblock0, vi);
    
    for j = setdiff(1:n, i)
        idxsj = binmap{j};
        vj = gauge{j}(:);
        vj = vj/sum(vj);
        
        [block, ~, rowavg, avg] = gaugify(params.couplings(idxsi, idxsj), vi, vj);
        newcouplings(idxsi, idxsj) = block;
        newcouplings(idxsj, idxsi) = block';
        
        dblock = dblock + 2*(rowavg - avg);
    end
    
    newcouplings(idxsi, idxsi) = diag(dblock);
end

newparams.couplings = newcouplings;

end

function [res, colavg, rowavg, avg] = gaugify(mat, vi, vj)
% Return a matrix in the specified gauge, as well as the vectors of average
% values per row and column weighted by vi/vj, and the average of all the
% elements in the block.

colavg = vi(:)'*mat;
rowavg = mat*vj(:);
avg = dot(colavg, vj);

res = mat - repmat(colavg, size(mat, 1), 1) - repmat(rowavg, 1, size(mat, 2)) + avg;

end