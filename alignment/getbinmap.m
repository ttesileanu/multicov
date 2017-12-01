function map = getbinmap(alignment)
% GETBINMAP Get map from alignment residues to groups of binary alignment
% columns.
%   map = GETBINMAP(alignment) gets a cell array of vectors. The vector
%   at map{i} is made up of the positions in the binary alignment
%   corresponding to residue index i in the alignment. The argument to the
%   function can be either an alignment or a binary alignment or, more
%   generally, a structure with fields alphabets and alphawidths.

% Tiberiu Tesileanu (2012-2014)

n = sum(alignment.alphawidths);
map = cell(n, 1);

nalpha = length(alignment.alphabets);
itot = 1;
ibin = 1;
for a = 1:nalpha
    width = alignment.alphawidths(a);
    alphabet = alignment.alphabets{a};
    nlett = length(alphagetletters(alphabet, 'nogap'));
    for i = 1:width
        map{itot} = ibin:(ibin + nlett - 1);
        ibin = ibin + nlett;
        itot = itot + 1;
    end
end

end