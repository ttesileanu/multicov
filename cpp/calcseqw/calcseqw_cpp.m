% CALCSEQW_CPP Internal function.
%
% Get adjacency matrix of the graph induced on an alignment by sequence
% similarity.
% (C++ helper function)
%   There is little reason to call this function directly. Instead, use
%   esimatseqw or eliminatesimilar, which automatically call the C++ function
%   when it is available.
%
%   adjmat = CALCSEQW_CPP(alignmatrix, max_seqid) returns a lower-triangular,
%   sparse adjacency matrix on the sequences in the given alignment matrix,
%   where an edge between two sequences means that they are at least max_seqid
%   identical in terms of normalized Hamming distance. The alignment should be
%   a character matrix.
%
%   See also: ESTIMATESEQW, ELIMINATESIMILAR.

% Tiberiu Tesileanu (2013-2014)
