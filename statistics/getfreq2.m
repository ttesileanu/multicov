function freq2 = getfreq2(alignment)
% GETFREQ2 Get the second moment of the distribution for the alignment.
%   freq2 = GETFREQ2(alignment) returns the second moment of the empirical
%   distribution implied by the alignment -- that is, the joint
%   frequencies of pairs of residues.
%
% See also GETFREQ1.

% Tiberiu Tesileanu (2012-2014)

if alncheck(alignment)
    binalign = alntobin(alignment);
elseif bincheck(alignment)
    binalign = alignment;
elseif statscheck(alignment)
    if isfield(alignment, 'freq2')
        freq2 = alignment.freq2;
    else
        freq2 = alignment.cmat + alignment.freq1(:)*alignment.freq1(:)';
    end
    return;
else
    error([mfilename ':badaln'], 'The argument should be an alignment or a statistics structure.');
end

freq2 = full(binalign.data'*diag(sparse(binalign.seqw(:)))*binalign.data/sum(binalign.seqw(:)));

end
