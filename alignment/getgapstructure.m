function gaps = getgapstructure(alignment)
% GETGAPSTRUCTURE Get a matrix that identifies the location of the gaps in
% the alignment.
%   gaps = GETGAPSTRUCTURE(alignment) returns a logical matrix with ones
%   in places where there are gaps in the alignment, zeros everywhere else.
%   For gapped alphabets, the structure returned is the same as if the
%   alphabets weren't gapped.

% Tiberiu Tesileanu (2012-2014)

if alncheck(alignment)
    ranges = getalpharanges(alignment);
    gaps = false(size(alignment.data));
    % do this one alphabet at a time
    for i = 1:size(ranges, 2)
        block = alignment.data(:, ranges(1, i):ranges(2, i));
        alphabet = alignment.alphabets{i};
        % need the ungapped version of the alphabets here
        if length(alphabet) > 3 && strcmp(alphabet(1:3), 'gap')
            alphabet = alphabet(4:end);
        end
        letters = alphagetletters(alphabet);
        % make sure to use the gap from the alphabet (assumed to be the first letter)
        gapletter = letters(1);
        gapblock = (block == gapletter);
        gaps(:, ranges(1, i):ranges(2, i)) = gapblock;
    end
elseif bincheck(alignment)
    binalign = alignment;
    
    Nseq = size(binalign.data, 1);
    Npos = sum(binalign.alphawidths);
    gaps = false(Nseq, Npos);

    alnmap = getbinmap(binalign);
    alphacumsum = cumsum(binalign.alphawidths);
    for i = 1:Npos
        block = binalign.data(:, alnmap{i});
        % gaps are represented by all zeros in the binary alignment
        % UNLESS this is a gapped alphabet!
        alphabet = binalign.alphabets{find(i <= alphacumsum, 1)};
        if length(alphabet) > 3 && strcmp(alphabet(1:3), 'gap')
            gapblock = block(:, 1);
        else
            gapblock = ~any(block, 2);
        end
        gaps(:, i) = gapblock;
    end
else
    error([mfilename ':badarg'], 'The first argument should an alignment structure.');
end

end
