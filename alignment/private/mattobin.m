function bin = mattobin(data, alphabet)
% MATTOBIN Internal function.
%
% Convert alignment to a binary alignment tensor.
%   bin = MATTOBIN(data, alphabet) converts the data matrix
%   into a binary alignment tensor, according to the given alphabet.
%   bintensor(k, a, i) = 1 iff sequence k has letter a at position i. In
%   this representation, gaps are not considered letters of the alphabet;
%   they are represented by bintensor(k, :, i) = 0. It then partially
%   flattens the tensor, returning a matrix in which the position and
%   character indices are combined.
%
%   See also ALNTOBIN.

% Tiberiu Tesileanu (2012, 2014)

letters = alphagetletters(alphabet, 'nogap');
nlett = length(letters);
bin = zeros([size(data, 1) nlett size(data, 2)]);
for i = 1:nlett
    bin(:, i, :) = (data == letters(i));
end
% reshape into a matrix
bin = reshape(bin, [size(bin, 1) size(bin, 2)*size(bin, 3)]);

end