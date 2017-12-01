function res = selectalpha(input, aidxs)
% SELECTALPHA Select some of the alphabets of an alignment or statistics
% structure.
%   res = SELECTALPHA(input, aidxs) generates a new alignment or statistics
%   structure that keeps only the alphabets indexed by 'aidxs'. Apart from
%   having elements from 1 to length(input.alphabets), there are no
%   restrictions on 'aidxs': alphabets can be repeated and their order can
%   be changed by using this function.

% Tiberiu Tesileanu (2014)

res = input;

% reordering the alphabets, alphawidths, and refseq information is easy
res.alphabets = input.alphabets(aidxs);
res.alphawidths = input.alphawidths(aidxs);
res.refseq = input.refseq(aidxs);

% working on the data is a bit more involved
if alncheck(input)
    ranges = getalpharanges(input);
    selranges = ranges(:, aidxs);
    idxs = cell2mat(arrayfun(@(i) flatten(selranges(1, i):selranges(2, i)), ...
        flatten(1:size(selranges, 2)), 'uniform', false));
    res.data = input.data(:, idxs);
else
    binranges = getalpharanges(input, 'binary');
    selbinranges = binranges(:, aidxs);
    binidxs = cell2mat(arrayfun(@(i) flatten(selbinranges(1, i):selbinranges(2, i)), ...
        flatten(1:size(selbinranges, 2)), 'uniform', false));
    if bincheck(input)
        res.data = input.data(:, binidxs);
    elseif statscheck(input)
        res.freq1 = input.freq1(binidxs);
        res.cmat = input.cmat(binidxs, binidxs);
        if isfield(input, 'freq2')
            res.freq2 = input.freq2(binidxs, binidxs);
        end
    else
        error([mfilename ':badinput'], 'The input must be an alignment or a statistics structure.');
    end
end

end