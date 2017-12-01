function newalign = alnscramble(alignment)
% ALNSCRAMBLE Scramble alignment data to remove correlations.
%   newalign = ALNSCRAMBLE(alignment) scrambles the alignment by shuffling
%   the letters in each column independently of any other columns. This
%   effectively removes correlations between positions without affecting
%   positional frequencies.
%
%   Note that this only works properly if the sequence weights are trivial
%   for the alignment. Otherwise, positional frequencies may change. The
%   function keeps the sequence weights vector of the original alignment.

% Tiberiu Tesileanu (2012-2014)

if ~alncheck(alignment)
    error([mfilename ':badarg'], 'First argument should be an alignment structure.');
end

[Nseq, Npos] = size(alignment.data);

% make sure we have all the tags for a character alignment
newalign = alnmake;
% intentionally don't copy extra parameters (such as 'annotations') which
% would be meaningless for the scrambled alignment
newalign.alphabets = alignment.alphabets;
newalign.alphawidths = alignment.alphawidths;
newalign.seqw = alignment.seqw;
newalign.data = alignment.data;
newalign.refseq = alignment.refseq;
for i = 1:Npos
    newalign.data(:, i) = alignment.data(randperm(Nseq), i);
end

end
