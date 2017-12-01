function [newalign, changed] = alnincludegaps(alignment)
% ALNINCLUDEGAPS Return a new alignment where gaps are treated identically
% to other letters.
%   newalign = ALNINCLUDEGAPS(alignment) returns a new alignment where gaps
%   are treated identically to other letters. This is essentially done by
%   replacing the alphabets by 'gapped' alphabets where the first letter
%   (usually considered the gap) is a dummy that does not appear in the
%   alignment data, while the gap appears as a normal letter.
%
%   newbinalign = ALNINCLUDEGAPS(binalign) also works on binary alignments.
%   Apart from changing to 'gapped' alphabets, this function also updates
%   the binary data by including columns for gaps.
%
%   Note that apart from binary alignments, this function treats all
%   structures having an 'alphabets' field equally. See includegaps.m for a
%   function that performs the other needed transformations for things like
%   statistics structures.
%
%   [..., changed] = ALNINCLUDEGAPS(...) returns 'changed' equal to true if
%   a change needed to be made, false otherwise.
%
% See also: ALPHAGETLETTERS, INCLUDEGAPS.

% Tiberiu Tesileanu (2012-2014)

% turn the alphabets into gapped alphabets
[newalign, changed] = alphaincludegaps(alignment);

if changed && bincheck(alignment)
    binmap = getbinmap(alignment);
    newbinmap = getbinmap(newalign);
    alphasum = cumsum(alignment.alphawidths);
    newalign.data = zeros(size(alignment.data, 1), newbinmap{end}(end));
    for i = 1:length(binmap)
        idxsnew = newbinmap{i};
        idxsold = binmap{i};
        a = find(i <= alphasum, 1);
        if strcmp(newalign.alphabets{a}, alignment.alphabets{a})
            % no transformation needed here
            newalign.data(:, idxsnew) = alignment.data(:, idxsold);
        else
            newalign.data(:, idxsnew(2:end)) = alignment.data(:, idxsold);
            newalign.data(:, idxsnew(1)) = ~any(alignment.data(:, idxsold), 2);
        end
    end
end

end