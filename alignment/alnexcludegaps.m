function [newalign, changed] = alnexcludegaps(alignment)
% ALNEXCLUDEGAPS Return a new alignment where gaps are treated normally
% (i.e., ignored in binary alignments, statistics structures, etc.)
%   newalign = ALNEXCLUDEGAPS(alignment) returns a new alignment where gap
%   values are meant to be absent from binary alignments, statistics
%   structures, etc. calculated from this alignment. This is essentially
%   done by replacing the 'gapped' alphabets by their ungapped counterparts
%   (i.e., removing the prefix 'gap' from alphabets that have it).
%
%   newbinalign = ALNINCLUDEGAPS(binalign) also works on binary alignments.
%   Apart from changing to 'ungapped' alphabets, this function also updates
%   the binary data by excluding the gap columns.
%
%   Note that apart from binary alignments, this function treats all
%   structures having an 'alphabets' field equally. See excludegaps.m for a
%   function that performs the other needed transformations for things like
%   statistics structures.
%
%   [..., changed] = ALNEXCLUDEGAPS(...) returns 'changed' equal to true if
%   a change needed to be made, false otherwise.
%
% See also: ALPHAGETLETTERS, EXCLUDEGAPS, ALNINCLUDEGAPS.

% Tiberiu Tesileanu (2014)

% turn the alphabets into ungapped alphabets
[newalign, changed] = alphaexcludegaps(alignment);

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
            newalign.data(:, idxsnew) = alignment.data(:, idxsold(2:end));
        end
    end
end

end