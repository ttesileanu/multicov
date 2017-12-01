function freq1 = getfreq1(alignment)
% GETFREQ1 Get the frequencies (first moment of the distribution) for
% an alignment.
%   freq1 = GETFREQ1(alignment) returns a vector of positional frequencies
%   of the residues in the alignment. The argument can be either a
%   character alignment, a binary alignment, or a statistics structure (in
%   which case freq1 is directly returned from the structure).
%
% See also GETFREQ2.

% Tiberiu Tesileanu (2012-2014)

if alncheck(alignment)
    binalign = alntobin(alignment);
elseif bincheck(alignment)
    binalign = alignment;
elseif statscheck(alignment)
	freq1 = alignment.freq1;
    return;
else
    error([mfilename ':badaln'], 'The argument should be an alignment or a statistics structure.');
end

freq1 = binalign.seqw(:)'*binalign.data/sum(binalign.seqw);
freq1 = freq1(:);

end