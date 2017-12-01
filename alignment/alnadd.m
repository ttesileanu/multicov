function align = alnadd(alignment, data, alphabet)
% ALNADD Place two alignments next to each other.
%   align = ALNADD(alignment, other_alignment) generates a new alignment
%   that contains the two alignments placed next two each other. An error
%   is generated if the number of sequences is not the same in both of
%   them. Note that the annotations and sequence weights from the first
%   alignment are kept (except in the degenerate case in which this
%   alignment is empty).
%
%   align = ALNADD(alignment, data, alphabet) adds the data (presented as a
%   matrix) to the right of the block of data already contained by the
%   alignment, and returns the new alignment that is formed in this way. An
%   error is generated if the number of sequences does not match. This is
%   equivalent to ALNADD(alignment, alnmake(data, alphabet)).
%
%   See also ALNMAKE, ALNCHECK.

% Tiberiu Tesileanu (2012-2014)

% check the input
if ~alncheck(alignment)
    error([mfilename ':badaln'], 'The first argument should be an alignment structure.');
end
is_data_aln = alncheck(data);
if ~ismatrix(data) && ~is_data_aln
    error([mfilename ':baddata'], 'The second argument should be a matrix or an alignment structure.');
end
if ~is_data_aln && nargin < 3
    error([mfilename ':noalpha'], 'When the second argument is a matrix, the third argument (the alphabet) is needed.');
end

N0 = size(alignment.data, 1);

% make sure all fields are copied (e.g., 'annotations' and 'seqw')
align = alignment;

% handle the first form of the command
if ~alncheck(data)
    other_align = alnmake(data, alphabet);
else
    other_align = data;
end

% concatenating two alignments
if N0 == 0
    % if the first alignment was empty, just return the alignment we're adding
    align = other_align;
elseif ~isempty(other_align.data)
    % if the alignment we're adding is empty, there's nothing to do
    if size(alignment.data, 1) ~= size(other_align.data, 1)
        error([mfilename ':badalnsize'], 'Mismatch between size of the two alignments.');
    end
    
    align.alphabets = [align.alphabets other_align.alphabets];
    align.alphawidths = [align.alphawidths other_align.alphawidths];
    align.data = [align.data other_align.data];
    
    align.refseq = [align.refseq(:) ; other_align.refseq(:)];
end

end