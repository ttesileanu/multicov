function newbinalign = binproject(binalign, dirvec)
% BINPROJECT Project each position in a binary alignment onto a vector.
%   newbinalign = BINPROJECT(binalign, dirvec) returns a binary alignment
%   in which the columns for each alignment position are projected onto a
%   corresponding vector from dirvec. If, for example, binalign is a
%   single-alphabet protein alignment, then columns 1:20 would be projected
%   onto dirvec(1:20), columns 21:40 onto dirvec(21:40), etc. The
%   projection is done by taking the dot product with the (L^2) normalized
%   vector.

% Tiberiu Tesileanu (2012-2014)

if ~bincheck(binalign)
    error([mfilename ':badaln'], 'The first argument should be a binary alignment.');
end

binmap = getbinmap(binalign);
newdata = zeros(size(binalign.data, 1), sum(binalign.alphawidths));
for i = 1:size(newdata, 2)
    idxs = binmap{i};
    subvec = dirvec(idxs);
    nvec = norm(subvec);
    if nvec > 0
        subvec = subvec / nvec;
    end
    newdata(:, i) = binalign.data(:, idxs)*subvec;
end

% copy all the data first; this keeps any data that the alignment might
% contain that we don't know about (and takes no time because Matlab
% implements lazy copying)
newbinalign = binalign;

newbinalign.data = newdata;
newbinalign.alphabets = {'binary'};
newbinalign.alphawidths = size(newbinalign.data, 2);

% XXX store some information about the original alignment?
% newbinalign.consensus = getconsensus(binalign);
% newbinalign.consalphabets = binalign.alphabets;
% newbinalign.consalphawidths = binalign.alphawidths;

end