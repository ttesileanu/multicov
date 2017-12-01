function isbin = bincheck(binalign)
% BINCHECK True if input is a binary alignment structure.
%   BINCHECK(binalign) returns true if binalign is a binary alignment
%   structure. This should contain fields 'type', 'alphabets', 'alphawidths',
%   'data', 'seqw', 'refseq', and 'annotations'. All except for 'data' and 
%   'type' have the same meanings as for a character aligment; see alnmake
%   for a description. 'Type' should be 'binary'; data is a numeric matrix,
%   normally a sparse matrix containing the binary alignment: each position
%   in the binary alignment corresponds to a certain letter X of the alphabet
%   at a certain position k in the ith sequence; it is nonzero only if the
%   alignment contains letter X at position k in sequence i. The nonzero
%   elements are 1 in the absence of positional weights.
%
%   See also ALNMAKE, ALNCHECK, STATSCHECK.

% Tiberiu Tesileanu (2012-2014)

isbin = (genalncheck(binalign) && strcmp(binalign.type, 'binary') && ...
    ismatrix(binalign.data) && ...
    (isnumeric(binalign.data) || islogical(binalign.data)));

end
