function newdist = catdiststruct(dist1, dist2)
% CATDISTSTRUCT Concatenate two distance structures.
%   newdist = CATDISTSTRUCT(dist1, dist2) makes a new distance structure
%   that contains dist1 and dist2 in that order. The distances between
%   structures in dist1 and structures in dist2 are marked as unknown by
%   setting the appropriate entries in the 'validmask' field to false.
%
% See also: GETDISTSTRUCT.

% Tiberiu Tesileanu (2014)

n1 = size(dist1.distmat, 1);
n2 = size(dist2.distmat, 1);
newdist.distmat = [dist1.distmat zeros(n1, n2) ; zeros(n2, n1) dist2.distmat];
newdist.validmask = [dist1.validmask false(n1, n2) ; false(n2, n1) dist2.validmask];
newdist.refseq = [dist1.refseq(:) ; dist2.refseq(:)];

end