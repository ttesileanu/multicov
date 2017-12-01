function res = selectcol(input, idxs)
% SELECTCOL Select some of the columns of an alignment or statistics structure.
%   res = SELECTCOL(input, idxs) selects ony the columns given in 'idxs'
%   from the (binary or character) alignment or statistics structure
%   'input'. The indices are given in terms of complete sequence residues
%   (i.e., index 1 would mean positions 1:20 for a protein alignment). The
%   indices can also be given as a logical array of length sum(input.alphawidths).
%
%   The 'idxs' vector should not have repeated values, and all the indices
%   corresponding to a given alignment should be adjacent to each other. It
%   is, however, possible to switch columns around within alphabets by using
%   this function. The results are undefined if trying to switch indices
%   from different alignments.
%
%   See also: SELECTSEQ.

% Tiberiu Tesileanu (2012-2014)

ranges = getalpharanges(input);

% make sure to copy any extra data that the alignment might contain
% (e.g., annotations)
res = input;

if islogical(idxs) || isequal(flatten(unique(idxs)), [0 ; 1])
    idxs = find(idxs);
end

if length(idxs) ~= length(unique(idxs))
    error([mfilename ':repidx'], 'Indices should not repeat.');
end

% updating the alphabets is more involved... make a mask of the
% positions we keep
mask = false(1, sum(input.alphawidths));
mask(idxs) = true;

newalphabets = {};
newalphawidths = [];
newrefseq = [];
crt = 1;
for i = 1:length(input.alphabets)
    submask = mask(ranges(1, i):ranges(2, i));
    newwidth = sum(submask);
    if newwidth > 0
        newalphabets = [newalphabets input.alphabets{i}]; %#ok<AGROW>
        newalphawidths = [newalphawidths newwidth]; %#ok<AGROW>
        crtrefseq = input.refseq(i);
        subidxs = idxs(idxs >= ranges(1, i) & idxs <= ranges(2, i)) - ranges(1, i) + 1;
        crtrefseq.map = input.refseq(i).map(subidxs);
        newrefseq = [newrefseq crtrefseq]; %#ok<AGROW>
        crt = crt + 1;
    end
end

if alncheck(input)
    % updating the data for character alignments is easy
    res.data = input.data(:, idxs);
else
    binmap = getbinmap(input);
    binidxs = cell2mat(flatten(cellfun(@flatten, binmap(idxs), 'uniform', false)));
    if bincheck(input)
        res.data = input.data(:, binidxs);
    elseif statscheck(input)
        res.freq1 = input.freq1(binidxs);
        res.cmat = input.cmat(binidxs, binidxs);
        if isfield(input, 'freq2')
            res.freq2 = input.freq2(binidxs, binidxs);
        end
    else
        error([mfilename ':badinput'], 'The input should be a character or binary alignment, or a statistics structure.');
    end
end

res.alphabets = newalphabets;
res.alphawidths = newalphawidths;
res.refseq = newrefseq;

end