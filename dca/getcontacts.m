function contacts = getcontacts(redJ, structure, varargin)
% GETCONTACTS Find contact predictions from DCA results.
%   contacts = GETCONTACTS(redJ, structure) uses the reduced coupling
%   matrix redJ to make contact predictions. The 'alphabets' and
%   'alphawidths' fields of 'structure' show how 'redJ' is split between
%   the different alphabets. The 'refseq' field identifies the mapping to
%   reference sequences. See the options for changing the number of
%   contacts predicted and how these are distributed.
%
%   The output is a structure with the following fields:
%    'mapped'
%       A cell matrix of size length(structure.alphabets). The cell at
%       position (i, j) identifies the contacts between residues in the ith
%       alphabet and residues in the jth alphabet. This is done using a
%       2 x k cell array in which the top row corresponds to the ith
%       alphabet and the bottom one corresponds to the jth alphabet. The
%       entries correspond to reference sequence positions, as identified
%       by refseq. The need for a cell array arises from the fact that the
%       reference sequence mapping might be done using strings and thus may
%       not be numeric.
%    'matrix'
%       A matrix of the same size as redJ showing where the predicted
%       contacts are. This is essentially a thresholded version of redJ.
%    'rawpairs'
%       A 2 x k vector identifying the indices in 'matrix' where the k
%       predicted contacts occur. The top value is always smaller than the
%       bottom one.
%    'redJ'
%       A copy of the matrix that was used for the prediction.
%    'refseq'
%       A structure analogous to that encountered in alignments and
%       statistics structures, identifying the mapping from indices in
%       'redJ' to some reference sequences. This is taken from 'structure'.
%
%   Options:
%    'seqcutoff' <n/v>
%       If non-zero, ignore contact pairs in which the residues are <= n
%       away. Note that the distance is measured in terms of refseq-mapped
%       indices when the refseq mapping is numeric, but in terms of raw
%       indices if it is not. If this is a scalar, the same threshold is
%       used for all alphabets. Otherwise, the ith threshold is used for
%       the ith alphabet.
%       (default: 4)
%    'n' <n/m>
%       This can either be a number or a square matrix with size equal to
%       the number of alphabets in 'structure'. When it is a scalar, the
%       largest n entries in redJ are assumed to be contacts, regardless of
%       which alphabets they are in. When it is a matrix m, m(i, j) gives
%       the number of contacts between alphabet i and alphabet j to find.
%       (default: equal to size of 'redJ')
%
% See also: GETDCA, GETMAXENTPARAMS, REDUCECOUPLINGS.

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('seqcutoff', 4, @(x) isnumeric(x) && isvector(x));
parser.addParamValue('n', [], @(x) isnumeric(x) && ismatrix(x) && size(x, 1) == size(x, 2));

% parse
parser.parse(varargin{:});
params = parser.Results;

if ~isnumeric(redJ) || ~ismatrix(redJ) || size(redJ, 1) ~= size(redJ, 2)
    error([mfilename ':badredJ'], 'redJ should be a square numeric matrix.');
end
if size(redJ, 1) ~= sum(structure.alphawidths)
    error([mfilename ':badstruc'], 'redJ size does not match structure.');
end

% set default n
if isempty(params.n)
    params.n = size(redJ, 1);
end

% handle scalar seqcutoff
nalphas = length(structure.alphabets);
if isscalar(params.seqcutoff)
    params.seqcutoff = repmat(params.seqcutoff, nalphas, 1);
end

% sort the elements in the matrix and find the valid contacts (the ones
% obeying the seqcutoff)
[~, sortorder] = sort(redJ(:), 'descend');
[sorti, sortj] = ind2sub(size(redJ), sortorder);
mask = (sorti < sortj);
sorti = sorti(mask);
sortj = sortj(mask); % this automatically gets rid of 'contacts' between a residue and itself
% this maps each matrix index to the corresponding alphabet
alphamap = cell2mat(arrayfun(@(i) i*ones(structure.alphawidths(i), 1), (1:length(structure.alphabets))', 'uniform', false));
sorti_alpha = alphamap(sorti);
sortj_alpha = alphamap(sortj);

% impose seqcutoff
mask = true(size(sorti));
alphacumsum = [1 ; 1 + flatten(cumsum(structure.alphawidths))];
for a = 1:nalphas
    % can only do it for numeric mappings
    if isnumeric(structure.refseq(a).map)
        samealpha_mask = (sorti_alpha == a & sortj_alpha == a);
        alphastart = alphacumsum(a);
        mapped_i = structure.refseq(a).map(sorti(samealpha_mask) - alphastart + 1);
        mapped_j = structure.refseq(a).map(sortj(samealpha_mask) - alphastart + 1);
        submask = (abs(mapped_j - mapped_i) > params.seqcutoff(a));
        mask(samealpha_mask) = submask;
    end
end

sorti = sorti(mask);
sortj = sortj(mask);

sorti_alpha = sorti_alpha(mask);
sortj_alpha = sortj_alpha(mask);

if isscalar(params.n)
    % impose the number cutoff independent of alphabet
    if length(sorti) > params.n
        sorti = sorti(1:params.n);
        sortj = sortj(1:params.n);

        sorti_alpha = sorti_alpha(1:params.n);
        sortj_alpha = sortj_alpha(1:params.n);
    elseif length(sorti) < params.n
        warning([mfilename ':missing'], ['Couldn''t find all the required ' int2str(params.n) ' contacts.']);
    end
end

contacts.mapped = cell(nalphas, nalphas);
mask = false(size(sorti));
for a = 1:nalphas
    alphai_a_mask = (sorti_alpha == a);
    alphastart_a = alphacumsum(a);
    for b = 1:nalphas
        alphaj_b_mask = (sortj_alpha == b);
        alphastart_b = alphacumsum(b);
        
        % a pair (i, j) implies the pair (j, i), but only one of them will
        % be present in [sorti, sortj]! we need to search for both
        alpha_idxs = find(...
            (alphai_a_mask & alphaj_b_mask) | (sorti_alpha == b & sortj_alpha == a));
        
        if ~isscalar(params.n)
            local_n = params.n(a, b);
            if length(alpha_idxs) > local_n
                alpha_idxs = alpha_idxs(1:local_n);
            elseif length(alpha_idxs) < local_n
                warning([mfilename ':missing2'], ['Couldn''t find all the required ' int2str(local_n) ' contacts for the alphabet pair ' int2str(a) ', ' int2str(b) '.']);
            end
        end
        
        subi = sorti(alpha_idxs);
        subj = sortj(alpha_idxs);

        % for the mapping, we need subi to contain indices in the a
        % alphabet, subj to contain indices in the b alphabet
        swap_mask = (~alphai_a_mask(alpha_idxs));
        
        subj_tmp = subj(swap_mask);
        subj(swap_mask) = subi(swap_mask);
        subi(swap_mask) = subj_tmp;
        
        mappedi = structure.refseq(a).map(subi - alphastart_a + 1);
        mappedj = structure.refseq(b).map(subj - alphastart_b + 1);
        
        if isnumeric(mappedi)
            mappedi = num2cell(mappedi);
        end
        if isnumeric(mappedj)
            mappedj = num2cell(mappedj);
        end
        
        contacts.mapped{a, b} = [mappedi(:)' ; mappedj(:)'];
        
        mask(alpha_idxs) = true;
    end
end

if ~isscalar(params.n)
    sorti = sorti(mask);
    sortj = sortj(mask);
end

contacts.rawpairs = [sorti(:)' ; sortj(:)'];

contacts.matrix = false(size(redJ));
contacts.matrix(sub2ind(size(contacts.matrix), contacts.rawpairs(1, :), contacts.rawpairs(2, :))) = true;
contacts.matrix(sub2ind(size(contacts.matrix), contacts.rawpairs(2, :), contacts.rawpairs(1, :))) = true;

contacts.redJ = redJ;
contacts.refseq = structure.refseq;

end