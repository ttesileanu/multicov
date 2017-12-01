function map = findseqmap(seq1, seq2, alphabet)
% FINDSEQMAP Find mapping between two sequences.
%   map = FINDSEQMAP(seq1, seq2, alphabet) outputs an array mapping
%   positions in the first sequence to positions in the second one. Map is
%   as long as the first sequence and each entry in map is either 0 (if no
%   match to the second sequence was found), or equal to the matching
%   position in the second sequence. The alphabet can be 'protein',
%   'gapprotein', or 'aa' for proteins, and 'dna'/'rna', 'gapdna'/'gaprna',
%   or 'nt' for nucleic acids. Other alphabets are not supported at this
%   time.
%
%   map = FINDSEQMAP(alignres) takes in the 'alignment' output argument of
%   nwalign instead of the sequences themselves.
%
%   The results of this function are undefined if either sequence contains
%   gaps.
%
% See also: NWALIGN.

% Tiberiu Tesileanu (2014)

if nargin > 1
    switch alphabet
        case {'protein', 'aa', 'gapprotein'}
            alphabet = 'aa';
        case {'dna', 'rna', 'nt', 'gapdna', 'gaprna'}
            alphabet = 'nt';
        otherwise
            error([mfilename ':badalpha'], 'Unsupported alphabet.');
    end
    [~, alignres] = nwalign(seq1, seq2, 'glocal', true, 'alphabet', alphabet);
else
    alignres = seq1;
end

% find where the letters are in the first sequence -- that's where we need
% maps for
letts1 = isletter(alignres(1, :));
% indices in the second sequence are given by the number of non-gaps up to
% a given point
letts2 = isletter(alignres(3, :));
idxs2 = cumsum(letts2);
idxs2(~letts2) = 0;
map = idxs2(letts1);

% make sure we output a column vector, for consistency
map = map(:);

end